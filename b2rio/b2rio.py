import pandas as pd
import numpy as np
import nibabel as nib
from neurolang.frontend import NeurolangPDL
import argparse

from nilearn import datasets, image


def run():
    parser = argparse.ArgumentParser()
    parser.add_argument("--brain_path", nargs=1, type=str, default=None)
    parser.add_argument("--n_folds", nargs='?', type=int, default=150)
    parser.add_argument("--resample", nargs='?', type=int, default=1)
    parser.add_argument("--frac_sample", nargs='?', type=int, default=0.7)
    parser.add_argument("--radius", nargs='?', type=int, default=4)
    parser.add_argument("--tfIdf", nargs='?', type=str, default='1e-3')
    parser.add_argument("--output_file", nargs=1, type=str, default=None)
    value = parser.parse_args()

    # %%
    if value.brain_path is None:
        print('You need to provide a nifti image using the --brain_path argument')
        return

    if value.output_file is None:
        print('You need to provide a name for the output file using the --output_file argument')
        return

    brain_path = value.brain_path[0]
    n_folds = value.n_folds
    resample = value.resample
    radius = value.radius
    frac_sample = value.frac_sample
    tfIdf = value.tfIdf
    output_file = value.output_file[0]


    print('Starting analysis with the following parameters:')
    print(f'  brain_path = {brain_path}')
    print(f'  n_folds = {n_folds}')
    print(f'  resample = {resample}')
    print(f'  radius = {radius}')
    print(f'  tfIdf = {tfIdf}')
    print(f'  frac_sample = {frac_sample}')
    print(f'  output_file = {output_file}')

    mni_t1 = nib.load(datasets.fetch_icbm152_2009()['t1'])
    mni_t1 = image.resample_img(mni_t1, np.eye(3) * resample)

    pmaps_4d = image.resample_img(
        image.load_img(brain_path), mni_t1.affine, interpolation='nearest'
    )

    brain_regions_data = []
    regions2analyse = set()
    brain_data = pmaps_4d.dataobj
    non_zero = np.nonzero(pmaps_4d.dataobj)
    if len(pmaps_4d.dataobj.shape) == 4:
        for x, y, z, r in zip(*non_zero):
            label = int(brain_data[x][y][z][r])
            regions2analyse.add(label)
            d = (x, y, z, label)
            brain_regions_data.append(d)
    elif len(pmaps_4d.dataobj.shape) == 3:
        regions2analyse.add(1)
        for x, y, z in zip(*non_zero):
            d = (x, y, z, 1)
            brain_regions_data.append(d)
    else:
        print('The nifti file must contain 3 or 4 dimensions')
        return


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
        brain_regions_data,
        name='brain_det'
    )

    with nl.scope as e:
        e.ontology_terms[e.onto_name, e.cp] = (
            hasTopConcept[e.uri, e.cp] &
            label[e.uri, e.onto_name]
        )

        e.lower_terms[e.term, e.cp] = (
            e.ontology_terms[e.onto_term, e.cp] &
            (e.term == word_lower[e.onto_term])
        )

        e.f_terms[e.d, e.t, e.cp] = (
            e.terms[e.d, e.t] &
            e.lower_terms[e.t, e.cp]
        )

        f_term = nl.query((e.d, e.t, e.cp), e.f_terms(e.d, e.t, e.cp))

    filtered = f_term.as_pandas_dataframe()[['d', 't']]
    nl.add_tuple_set(filtered.values, name='filtered_terms')

    if len(pmaps_4d.dataobj.shape) == 4:
        print('Starting analysis for regions:', regions2analyse)
    else:
        print('Starting analysis')

    for region in regions2analyse:
        try:
            with nl.scope as e:
                e.jbd[e.x, e.y, e.z, e.d] = (
                    e.activations[e.d, e.x, e.y, e.z] &
                    e.brain_det[e.x1, e.y1, e.z1, e.region] &
                    (e.dist == e.EUCLIDEAN(e.x, e.y, e.z, e.x1, e.y1, e.z1)) &
                    (e.dist < radius) &
                    (e.region == region)
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
            return

        df = res.as_pandas_dataframe()

        df_cp = f_term.as_pandas_dataframe()[['cp', 't']]
        df_cp = df_cp.drop_duplicates()

        df = df.set_index('term').join(df_cp.set_index('t'))
        df.reset_index(inplace=True)
        df['cp'] = df.cp.fillna('')
        df = df.rename(columns={'cp': 'topConcept', 'index': 'term'})

        if len(regions2analyse) > 1:
            df.to_csv(f'{output_file}_region{region}.csv', index=False)
        else:
            df.to_csv(f'{output_file}.csv', index=False)

        print(f'Results ready!')

def run_probabilistic():
    parser = argparse.ArgumentParser()
    parser.add_argument("--brain_path", nargs=1, type=str, default=None)
    parser.add_argument("--n_folds", nargs='?', type=int, default=150)
    parser.add_argument("--resample", nargs='?', type=int, default=1)
    parser.add_argument("--frac_sample", nargs='?', type=int, default=0.7)
    parser.add_argument("--radius", nargs='?', type=int, default=4)
    parser.add_argument("--tfIdf", nargs='?', type=str, default='1e-3')
    parser.add_argument("--output_file", nargs=1, type=str, default=None)
    value = parser.parse_args()

    # %%
    if value.brain_path is None:
        print('You need to provide a nifti image using the --brain_path argument')
        return

    if value.output_file is None:
        print('You need to provide a name for the uotput fiile using the --output_file argument')
        return

    brain_path = value.brain_path
    n_folds = value.n_folds
    resample = value.resample
    radius = value.radius
    frac_sample = value.frac_sample
    tfIdf = value.tfIdf
    output_file = value.output_file


    print('Starting analysis with the following parameters:')
    print(f'  n_folds = {n_folds}')
    print(f'  resample = {resample}')
    print(f'  radius = {radius}')
    print(f'  tfIdf = {tfIdf}')
    print(f'  frac_sample = {frac_sample}')
    print(f'  folder_results = {output_file}')

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

    brain_regions_data = []
    regions2analyse = set()
    brain_data = pmaps_4d.dataobj
    non_zero = np.nonzero(pmaps_4d.dataobj)
    if len(pmaps_4d.dataobj.shape) == 4:
        for x, y, z, r in zip(*non_zero):
            p = int(brain_data[x][y][z][r])
            regions2analyse.add(p)
            d = (p, x, y, z, r)
            brain_regions_data.append(d)
    elif len(pmaps_4d.dataobj.shape) == 3:
        regions2analyse.add(1)
        for x, y, z in zip(*non_zero):
            p = int(brain_data[x][y][z])
            d = (p, x, y, z, 1)
            brain_regions_data.append(d)
    else:
        print('The nifti file must contain 3 or 4 dimensions')
        return

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
        name='brain_det'
    )

    with nl.scope as e:
        e.ontology_terms[e.onto_name, e.cp] = (
            hasTopConcept[e.uri, e.cp] &
            label[e.uri, e.onto_name]
        )

        e.lower_terms[e.term, e.cp] = (
            e.ontology_terms[e.onto_term, e.cp] &
            (e.term == word_lower[e.onto_term])
        )

        e.f_terms[e.d, e.t, e.cp] = (
            e.terms[e.d, e.t] &
            e.lower_terms[e.t, e.cp]
        )

        f_term = nl.query((e.d, e.t, e.cp), e.f_terms(e.d, e.t, e.cp))

    filtered = f_term.as_pandas_dataframe()[['d', 't']]
    nl.add_tuple_set(filtered.values, name='filtered_terms')

    if len(pmaps_4d.dataobj.shape) == 4:
        print('Starting analysis for regions:', regions2analyse)
    else:
        print('Starting analysis')

    for region in regions2analyse:

        try:
            with nl.scope as e:
                (e.jbd @ e.p)[e.x, e.y, e.z, e.d] = (
                    e.activations[e.d, e.x, e.y, e.z] &
                    e.brain_det[e.p, e.x1, e.y1, e.z1, e.region] &
                    (e.dist == e.EUCLIDEAN(e.x, e.y, e.z, e.x1, e.y1, e.z1)) &
                    (e.dist < radius) &
                    (e.region == region)
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
            return

        df = res.as_pandas_dataframe()

        df_cp = f_term.as_pandas_dataframe()[['cp', 't']]
        df_cp = df_cp.drop_duplicates()

        df = df.set_index('term').join(df_cp.set_index('t'))
        df.reset_index(inplace=True)
        df['cp'] = df.cp.fillna('')
        df = df.rename(columns={'cp': 'topConcept', 'index': 'term'})

        if len(regions2analyse) > 1:
            df.to_csv(f'{output_file}_region{region}.csv')
        else:
            df.to_csv(f'{output_file}.csv')

    print(f'Results ready!')

    return f'Results ready at {output_file}.csv'