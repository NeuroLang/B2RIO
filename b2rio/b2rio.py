import pandas as pd
import numpy as np
import nibabel as nib
from neurolang.frontend import NeurolangPDL
import argparse

from nilearn import datasets, image


def run():

    print('Parsing args')
    parser = argparse.ArgumentParser()
    parser.add_argument("--brain_path", nargs='?', type=str, default=None)
    parser.add_argument("--n_folds", nargs='?', type=int, default=150)
    parser.add_argument("--resample", nargs='?', type=int, default=1)
    parser.add_argument("--frac_sample", nargs='?', type=int, default=0.7)
    parser.add_argument("--radius", nargs='?', type=int, default=4)
    parser.add_argument("--tfIdf", nargs='?', type=str, default='1e-3')
    parser.add_argument("--folder_results", nargs='?', type=str, default='./')
    value = parser.parse_args()

    # %%
    if value.brain_path is None:
        print('You need to provide a nifti image using the --brain_path argument')
        return

    brain_path = value.brain_path
    n_folds = value.n_folds
    resample = value.resample
    radius = value.radius
    frac_sample = value.frac_sample
    tfIdf = value.tfIdf
    folder_results = value.folder_results


    print('Starting analysis with the following parameters:')
    print(f'  n_folds = {n_folds}')
    print(f'  resample = {resample}')
    print(f'  radius = {radius}')
    print(f'  tfIdf = {tfIdf}')
    print(f'  frac_sample = {frac_sample}')

    print(f'  folder_results = {folder_results}')

    mni_t1 = nib.load(datasets.fetch_icbm152_2009()['t1'])
    mni_t1 = image.resample_img(mni_t1, np.eye(3) * resample)

    pmaps_4d = image.resample_img(
        image.load_img(brain_path), mni_t1.affine, interpolation='nearest'
    )

    brain_regions_prob = []
    prob_region_data = pmaps_4d.dataobj
    non_zero = np.nonzero(pmaps_4d.dataobj)
    for x, y, z, r in zip(*non_zero):
        p = prob_region_data[x][y][z][r]
        d = (p, x, y, z)
        brain_regions_prob.append(d)

    ns_database_fn, ns_features_fn = datasets.utils._fetch_files(
        datasets.utils._get_dataset_dir('neurosynth'),
        [
            (
                'database.txt',
                'https://github.com/neurosynth/neurosynth-data/raw/master/current_data.tar.gz',
                {'uncompress': True}
            ),
            (
                'features.txt',
                'https://github.com/neurosynth/neurosynth-data/raw/master/current_data.tar.gz',
                {'uncompress': True}
            ),
        ]
    )

    ns_database = pd.read_csv(ns_database_fn, sep=f'\t')
    ijk_positions = (
        nib.affines.apply_affine(
            np.linalg.inv(mni_t1.affine),
            ns_database[['x', 'y', 'z']]
        ).astype(int)
    )
    ns_database['i'] = ijk_positions[:, 0]
    ns_database['j'] = ijk_positions[:, 1]
    ns_database['k'] = ijk_positions[:, 2]

    ns_features = pd.read_csv(ns_features_fn, sep=f'\t')
    ns_terms = (
        pd.melt(
                ns_features,
                var_name='term', id_vars='pmid', value_name='TfIdf'
        )
        .query(f'TfIdf > {tfIdf}')[['pmid', 'term']]
    )
    ns_docs = ns_features[['pmid']].drop_duplicates()

    cogAt = datasets.utils._fetch_files(
        datasets.utils._get_dataset_dir('CogAt'),
        [
            (
                'cogat.xml',
                'http://data.bioontology.org/ontologies/COGAT/download?'
                'apikey=8b5b7825-538d-40e0-9e9e-5ab9274a9aeb&download_format=rdf',
                {'move': 'cogat.xml'}
            )
        ]
    )[0]

    nl = NeurolangPDL()
    nl.load_ontology(cogAt)

    subclass_of = nl.new_symbol(name='neurolang:subClassOf')
    label = nl.new_symbol(name='neurolang:label')
    hasTopConcept = nl.new_symbol(name='neurolang:hasTopConcept')

    @nl.add_symbol
    def word_lower(name: str) -> str:
        return name.lower()

    ns_doc_folds = pd.concat(
        ns_docs.sample(frac=frac_sample, random_state=i).assign(fold=[i] * (int((len(ns_docs)* frac_sample))+1))
        for i in range(n_folds)
    )
    doc_folds = nl.add_tuple_set(ns_doc_folds, name='doc_folds')

    ns_data = ns_database[ns_database.space == 'MNI'][['id', 'i', 'j', 'k']].values
    activations = nl.add_tuple_set(ns_data, name='activations')
    terms = nl.add_tuple_set(ns_terms.values, name='terms')
    docs = nl.add_uniform_probabilistic_choice_over_set(
            ns_docs.values, name='docs'
    )

    nl.add_tuple_set(
        brain_regions_prob,
        name='julich_brain_det'
    )

    with nl.scope as e:
        e.ontology_terms[e.onto_name] = (
            hasTopConcept[e.uri, e.cp] &
            label[e.uri, e.onto_name]
        )

        e.lower_terms[e.term] = (
            e.ontology_terms[e.onto_term] &
            (e.term == word_lower[e.onto_term])
        )

        e.f_terms[e.d, e.t] = (
            e.terms[e.d, e.t] &
            e.lower_terms[e.t]
        )

        f_term = nl.query((e.d, e.t), e.f_terms(e.d, e.t))

    filtered = f_term.as_pandas_dataframe()
    nl.add_tuple_set(filtered.values, name='filtered_terms')

    try:
        with nl.scope as e:
            e.jbd[e.x, e.y, e.z, e.d] = (
                e.activations[e.d, e.x, e.y, e.z] &
                e.julich_brain_det[e.x1, e.y1, e.z1] &
                (e.dist == e.EUCLIDEAN(e.x, e.y, e.z, e.x1, e.y1, e.z1)) &
                (e.dist < radius)
            )

            e.img_studies[e.d] = e.jbd[..., ..., ..., e.d]

            e.img_left_studies[e.d] = e.docs[e.d] & ~(e.img_studies[e.d])

            e.term_prob[e.t, e.fold, e.PROB[e.t, e.fold]] = (
                (
                    e.filtered_terms[e.d, e.t]
                ) // (
                    e.img_studies[e.d] &
                    e.doc_folds[e.d, e.fold] &
                    e.docs[e.d]
                )
            )

            e.no_term_prob[e.t, e.fold, e.PROB[e.t, e.fold]] = (
                (
                    e.filtered_terms[e.d, e.t]
                ) // (
                    e.img_left_studies[e.d] &
                    e.doc_folds[e.d, e.fold] &
                    e.docs[e.d]
                )
            )

            e.ans[e.term, e.fold, e.bf] = (
                e.term_prob[e.term, e.fold, e.p] &
                e.no_term_prob[e.term, e.fold, e.pn] &
                (e.bf == (e.p / e.pn))
            )

            res = nl.query((e.term, e.fold, e.bf), e.ans[e.term, e.fold, e.bf])
    except Exception as e:
        print(f'ERROR! : {e}')

    df = res.as_pandas_dataframe()
    df.to_csv(f'{folder_results}b2rio_results.csv')

    return f'Results ready at {folder_results}b2rio_results.csv'

def run_probabilistic():
    print('Parsing args')
    parser = argparse.ArgumentParser()
    parser.add_argument("--brain_path", nargs='?', type=str, default=None)
    parser.add_argument("--n_folds", nargs='?', type=int, default=150)
    parser.add_argument("--resample", nargs='?', type=int, default=1)
    parser.add_argument("--frac_sample", nargs='?', type=int, default=0.7)
    parser.add_argument("--radius", nargs='?', type=int, default=4)
    parser.add_argument("--tfIdf", nargs='?', type=str, default='1e-3')
    parser.add_argument("--folder_results", nargs='?', type=str, default='./')
    value = parser.parse_args()

    # %%
    if value.brain_path is None:
        print('You need to provide a nifti image using the --brain_path argument')
        return

    brain_path = value.brain_path
    n_folds = value.n_folds
    resample = value.resample
    radius = value.radius
    frac_sample = value.frac_sample
    tfIdf = value.tfIdf
    folder_results = value.folder_results


    print('Starting analysis with the following parameters:')
    print(f'  n_folds = {n_folds}')
    print(f'  resample = {resample}')
    print(f'  radius = {radius}')
    print(f'  frac_sample = {frac_sample}')
    print(f'  tfIdf = {tfIdf}')
    print(f'  folder_results = {folder_results}')

    mni_t1 = nib.load(datasets.fetch_icbm152_2009()['t1'])
    mni_t1 = image.resample_img(mni_t1, np.eye(3) * resample)

    pmaps_4d = image.resample_img(
        image.load_img(brain_path), mni_t1.affine, interpolation='nearest'
    )

    brain_regions_prob = []
    non_zero = np.nonzero(pmaps_4d.dataobj)
    for x, y, z, p in zip(*non_zero):
        d = (p, x, y, z)
        brain_regions_prob.append(d)

    ns_database_fn, ns_features_fn = datasets.utils._fetch_files(
        datasets.utils._get_dataset_dir('neurosynth'),
        [
            (
                'database.txt',
                'https://github.com/neurosynth/neurosynth-data/raw/master/current_data.tar.gz',
                {'uncompress': True}
            ),
            (
                'features.txt',
                'https://github.com/neurosynth/neurosynth-data/raw/master/current_data.tar.gz',
                {'uncompress': True}
            ),
        ]
    )

    ns_database = pd.read_csv(ns_database_fn, sep=f'\t')
    ijk_positions = (
        nib.affines.apply_affine(
            np.linalg.inv(mni_t1.affine),
            ns_database[['x', 'y', 'z']]
        ).astype(int)
    )
    ns_database['i'] = ijk_positions[:, 0]
    ns_database['j'] = ijk_positions[:, 1]
    ns_database['k'] = ijk_positions[:, 2]

    ns_features = pd.read_csv(ns_features_fn, sep=f'\t')
    ns_terms = (
        pd.melt(
                ns_features,
                var_name='term', id_vars='pmid', value_name='TfIdf'
        )
        .query(f'TfIdf > {tfIdf}')[['pmid', 'term']]
    )
    ns_docs = ns_features[['pmid']].drop_duplicates()

    cogAt = datasets.utils._fetch_files(
        datasets.utils._get_dataset_dir('CogAt'),
        [
            (
                'cogat.xml',
                'http://data.bioontology.org/ontologies/COGAT/download?'
                'apikey=8b5b7825-538d-40e0-9e9e-5ab9274a9aeb&download_format=rdf',
                {'move': 'cogat.xml'}
            )
        ]
    )[0]

    nl = NeurolangPDL()
    nl.load_ontology(cogAt)

    subclass_of = nl.new_symbol(name='neurolang:subClassOf')
    label = nl.new_symbol(name='neurolang:label')
    hasTopConcept = nl.new_symbol(name='neurolang:hasTopConcept')

    @nl.add_symbol
    def word_lower(name: str) -> str:
        return name.lower()

    ns_doc_folds = pd.concat(
        ns_docs.sample(frac=frac_sample, random_state=i).assign(fold=[i] * (int((len(ns_docs)* frac_sample))+1))
        for i in range(n_folds)
    )
    doc_folds = nl.add_tuple_set(ns_doc_folds, name='doc_folds')

    ns_data = ns_database[ns_database.space == 'MNI'][['id', 'i', 'j', 'k']].values
    activations = nl.add_tuple_set(ns_data, name='activations')
    terms = nl.add_tuple_set(ns_terms.values, name='terms')
    docs = nl.add_uniform_probabilistic_choice_over_set(
            ns_docs.values, name='docs'
    )

    nl.add_tuple_set(
        brain_regions_prob,
        name='julich_brain_det'
    )

    with nl.scope as e:
        e.ontology_terms[e.onto_name] = (
            hasTopConcept[e.uri, e.cp] &
            label[e.uri, e.onto_name]
        )

        e.lower_terms[e.term] = (
            e.ontology_terms[e.onto_term] &
            (e.term == word_lower[e.onto_term])
        )

        e.f_terms[e.d, e.t] = (
            e.terms[e.d, e.t] &
            e.lower_terms[e.t]
        )

        f_term = nl.query((e.d, e.t), e.f_terms(e.d, e.t))

    filtered = f_term.as_pandas_dataframe()
    nl.add_tuple_set(filtered.values, name='filtered_terms')

    try:
        with nl.scope as e:
            (e.jbd @ e.p)[e.x, e.y, e.z, e.d] = (
                e.activations[e.d, e.x, e.y, e.z] &
                e.julich_brain_det[e.p, e.x1, e.y1, e.z1] &
                (e.dist == e.EUCLIDEAN(e.x, e.y, e.z, e.x1, e.y1, e.z1)) &
                (e.dist < radius)
            )

            e.img_studies[e.d] = e.jbd[..., ..., ..., e.d]

            e.img_left_studies[e.d] = e.docs[e.d] & ~(e.img_studies[e.d])

            e.term_prob[e.t, e.fold, e.PROB[e.t, e.fold]] = (
                (
                    e.filtered_terms[e.d, e.t]
                ) // (
                    e.img_studies[e.d] &
                    e.doc_folds[e.d, e.fold] &
                    e.docs[e.d]
                )
            )

            e.no_term_prob[e.t, e.fold, e.PROB[e.t, e.fold]] = (
                (
                    e.filtered_terms[e.d, e.t]
                ) // (
                    e.img_left_studies[e.d] &
                    e.doc_folds[e.d, e.fold] &
                    e.docs[e.d]
                )
            )

            e.ans[e.term, e.fold, e.bf] = (
                e.term_prob[e.term, e.fold, e.p] &
                e.no_term_prob[e.term, e.fold, e.pn] &
                (e.bf == (e.p / e.pn))
            )

            res = nl.query((e.term, e.fold, e.bf), e.ans[e.term, e.fold, e.bf])
    except Exception as e:
        print(f'ERROR! : {e}')

    df = res.as_pandas_dataframe()
    df.to_csv(f'{folder_results}b2rio_results.csv')

    return f'Results ready at {folder_results}b2rio_results.csv'