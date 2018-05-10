import pickle
import time
import os

from indra.sources.isi.api import process_preprocessed
from indra.sources.isi.preprocessor import IsiPreprocessor


def abstracts_runtime():
    pfile = '/Users/daniel/Downloads/text_content_sample.pkl'
    dump = pickle.load(open(pfile, 'rb'))
    all_abstracts = dump['pubmed']

    # abstract_counts = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    # abstract_counts = [20, 40, 60, 80, 100]
    # abstract_counts = [20, 50, 100, 1000, None]
    # abstract_counts = [20, 50, 100, 200, 300, 400, 500]
    abstract_counts = [10]
    times = []
    for count in abstract_counts:
        print('==============================================================')
        print('Count:', count)
        print('==============================================================')
        with open('isi_experiment_log.txt', 'a') as f:
            f.write('Reading and processing ' + str(count) + 'abstracts\n')
        start_time = time.time()

        if count is None:
            abstract_subset = all_abstracts
        else:
            abstract_subset = all_abstracts[:count]
            assert(len(abstract_subset) == count)

        preprocessed_dir = 'preprocessed_' + str(count)
        os.mkdir(preprocessed_dir)
        preprocessor = IsiPreprocessor(preprocessed_dir)
        preprocessor.preprocess_abstract_list(abstract_subset)

        output_dir = 'output_' + str(count)
        os.mkdir(output_dir)
        ip = process_preprocessed(preprocessor, 8, output_dir)
        print('Statements:', ip.statements)

        output_pickle = 'isi_processed_abstracts_' + str(count) + '.pkl'
        pickle.dump(ip, open(output_pickle, 'wb'))

        # How long did that take?
        ellapsed_s = time.time() - start_time
        ellapsed_min = ellapsed_s / 60.0
        times.append(ellapsed_min)

        with open('isi_experiment_log.txt', 'a') as f:
            f.write('Times so far: ')
            f.write(repr(times))
            f.write('\n')
    print('Times:')
    print(times)


if __name__ == '__main__':
    abstracts_runtime()
