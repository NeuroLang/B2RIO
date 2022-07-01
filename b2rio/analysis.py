import pandas as pd
import statistics
from neurolang.frontend import NeurolangPDL
from nilearn import datasets

class Analyzer():

    def run(self, regions_path, results_path='./'):

        lines = []
        with open(regions_path) as f:
            lines = f.readlines()

        count = 0
        res = []
        for line in lines:
            if count == 0:
                count += 1
                continue
            splited = line.split(' ')
            res.append((splited[0], ' '.join(splited[1:-1])[1:], splited[-1][:-2]))

        regions = pd.DataFrame(res, columns=['r_number', 'r_name', 'hemis'])
        regions = regions.astype({'r_number': 'int64'})

        res2 = {}
        for id_region in regions.r_number.values:

            term = pd.read_hdf(f'{results_path}neuro_paper_ri_term_probs_region{id_region}.hdf', key=f'results')

            res_t = []
            for t in term.t.unique():
                t1 = term[term.t == t]

                bayes_mode = statistics.mode(t1.bf.values)
                bayes_std = t1.bf.values.std()

                res_t.append((t, bayes_mode, bayes_std))



            res_reg = pd.DataFrame(res_t, columns=['term', 'bayes mode', 'bayes std'])#,  'log bayes mode', 'log bayes std'])
            res2[id_region] = res_reg


        df_res = pd.DataFrame([], columns=['term', 'bayes mode', 'bayes std', 'region', 'hemis'])
        for id_region, name, hemis in regions.values:
            df_region = res2[id_region].sort_values(['bayes mode'], ascending=False)
            df_region['region'] = name
            df_region['hemis'] = hemis


            df_res = df_res.append(df_region)

        cogAt_old = datasets.utils._fetch_files(
            datasets.utils._get_dataset_dir('ontology'),
            [
                (
                    'cogat_old.xml',
                    'https://data.bioontology.org/ontologies/COGAT/submissions/7/download?'
                    'apikey=8b5b7825-538d-40e0-9e9e-5ab9274a9aeb',
                    {'move': 'cogat_old.xml'}
                )
            ]
        )[0]

        nl = NeurolangPDL()
        nl.load_ontology(cogAt_old)

        part_of = nl.new_symbol(name='ro.owl:part_of')
        subclass_of = nl.new_symbol(name='neurolang:subClassOf')
        label = nl.new_symbol(name='neurolang:label')
        hasTopConcept = nl.new_symbol(name='neurolang:hasTopConcept')

        @nl.add_symbol
        def word_lower(name: str) -> str:
            return name.lower()

        # NeuroSynth
        _, ns_features_fn = datasets.utils._fetch_files(
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

        ns_features = pd.read_csv(ns_features_fn, sep=f'\t')
        ns_terms = (
            pd.melt(
                    ns_features,
                    var_name='term', id_vars='pmid', value_name='TfIdf'
            )
            .query('TfIdf > 1e-3')[['pmid', 'term']]
        )

        terms = nl.add_tuple_set(ns_terms.values, name='terms')

        with nl.scope as e:

            e.ontology_terms[e.cp, e.onto_name] = (
                hasTopConcept[e.uri, e.cp] &
                label[e.uri, e.onto_name]
            )

            e.lower_terms[e.cp, e.term] = (
                e.ontology_terms[e.cp, e.onto_term] &
                (e.term == word_lower[e.onto_term])
            )

            e.filtered_terms[e.cp, e.t] = (
                e.terms[..., e.t] &
                e.lower_terms[e.cp, e.t]
            )

            f_term = nl.query((e.topConcept, e.t), e.filtered_terms[e.topConcept, e.t])

        df_cp = f_term.as_pandas_dataframe()[['topConcept', 't']]
        df_cp = df_cp.drop_duplicates()
        df_res2 = df_res.set_index('term').join(df_cp.set_index('t'))
        df_res2 = df_res2.rename(columns={'cp': 'topConcept'})

        df_res2.to_hdf(f'{results_path}B2RIO results.hdf', key='results')
