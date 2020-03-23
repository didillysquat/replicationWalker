import requests
import pandas as pd
from bs4 import BeautifulSoup
import os
import sys
import compress_pickle
from multiprocessing import Pool, current_process


# TODO in this version of the class we will esentially use the current class as written passing in different
# base directories and using mp friendly vairbales to collect the data

class ReplicationWalkerHandler:
    def __init__(self, marker):
        self.marker = marker
        if self.marker == 'its2':
            self.remote_base_dir = "https://www.genoscope.cns.fr/sadc/tarapacific/METABARCODING/ITS2/ITS2_SYM_VAR_5.8S2_SYM_VAR_REV/"
        elif self.marker == '18s':
            self.remote_base_dir = "https://www.genoscope.cns.fr/sadc/tarapacific/METABARCODING/18S_V9/18S_V9_1389F_1510R/"
        elif self.marker == '16s_45':
            self.remote_base_dir = "https://www.genoscope.cns.fr/sadc/tarapacific/METABARCODING/16S_V4V5/Fuhrman_primers/"
        elif self.marker == '16s_full_45':
            self.remote_base_dir = "https://www.genoscope.cns.fr/sadc/tarapacific/METABARCODING/16S_Full_Length_plus_16S_V4V5/16S_FL_27F_1492R_plus_Fuhrman_primers/"
        self.exe_path = os.path.dirname(os.path.abspath(sys.argv[0]))
        self.authorisation_tup = self._make_auth_tup()
        self.headers = {'User-Agent': 'Benjamin Hume', 'From': 'benjamin.hume@kaust.edu.sa'}
        self.error_df_list_of_lists = None
        self.error_df_columns = ['barcode_id', 'readset_name', 'fwd_fastq_name', 'rev_fastq_name', 'pcr_sample_name', 'pcr_fl_sample_name' ,'dna_sample_name',
                                 'difference_category', 'color_code', 'directory']
        self.output_dir = "/home/humebc/projects/tara/replication_testing/output"
        os.makedirs(self.output_dir, exist_ok=True)

    def start_rep_walker_mp(self):
        # Get a list of the links that we will want to be visiting from the directory
        # below the remote_base_dir. These will be come the remote_base_dirs of the worker
        soup = BeautifulSoup(requests.get(self.remote_base_dir, auth=self.authorisation_tup, headers=self.headers).text, features="html.parser",)
        worker_base_dirs = [link.string for link in soup.find_all('a') if ((link.string not in ['Name', 'Last modified', 'Size', 'Description', 'Parent Directory', 'NEGATIVE_CONTROLS/', 'NEGATIVE_CONTROL/']) and ('/'.join([self.remote_base_dir.strip('/'), link.string]) not in self.remote_base_dir))]
        worker_base_dirs = [os.path.join(self.remote_base_dir, _) for _ in worker_base_dirs]
        with Pool(20) as p:
            self.error_df_list_of_lists = p.map(ReplicationWalkerWorker, worker_base_dirs)
        self.error_df_list_of_lists = [_.error_df_lists for _ in self.error_df_list_of_lists]
        # Now its time to create and write out our df
        data = []
        for _ in self.error_df_list_of_lists:
            data.extend(_)
        df = pd.DataFrame(data=data, columns=self.error_df_columns)
        df.to_csv(os.path.join(self.output_dir, f'replication_df_{self.marker}.csv'), index=False)

    def _make_auth_tup(self):
        auth_path = os.path.join(self.exe_path, 'auth.txt')
        with open(auth_path, 'r') as f:
            auth_lines = [line.rstrip() for line in f]
        return (auth_lines[0], auth_lines[1])


class ReplicationWalkerWorker:
    #TODO we need to modify this so that it takes in a single directory
    def __init__(self, remote_base_dir):
        self.remote_base_dir = remote_base_dir
        self.readset_info_dir = "/home/humebc/projects/tara/replication_testing/readset_csvs"
        self.readset_df = self._make_readset_info_dir()
        self.current_remote_dir = self.remote_base_dir
        self.error_df_lists = []
        self.done_list = set()
        self.done_and_empty_list = set()
        self.headers = {'User-Agent': 'Benjamin Hume', 'From': 'benjamin.hume@kaust.edu.sa'}
        self.exe_path = os.path.dirname(os.path.abspath(sys.argv[0]))
        self.authorisation_tup = self._make_auth_tup()
        self.fastq_gz_list_current = []
        self.links_to_visit_current = []
        self.last_fork_location = None
        self.home_dir_reached = None
        self.s = requests.Session()
        self.s.auth = self.authorisation_tup
        self.s.headers = self.headers
        self._walk()

    def _walk(self):
        while True:
            print(f'{current_process()}: Current directory: {self.current_remote_dir}')
            soup = BeautifulSoup(self.s.get(self.current_remote_dir).text, features="html.parser")
            self.links_to_visit_current = [link.string for link in soup.find_all('a') if
                              ((link.string not in ['Name', 'Last modified', 'Size', 'Description', 'Parent Directory', 'NEGATIVE_CONTROLS/']) and
                               ('/'.join([self.current_remote_dir.strip('/'), link.string]) not in self.done_list))]
            self.fastq_gz_list_current = [link.string for link in self.links_to_visit_current if 'fastq.gz' in link.string]
            self.links_to_visit_current = [link for link in self.links_to_visit_current if link not in self.fastq_gz_list_current]
            if len(self.links_to_visit_current) > 1:
                self.last_fork_location = self.current_remote_dir
            else:
                self.done_and_empty_list.add(self.current_remote_dir)
            if self.current_remote_dir == self.remote_base_dir and not self.links_to_visit_current:
                break

            if self.fastq_gz_list_current:
                try:
                    self._check_for_replicates()
                except IndexError as e:
                    print(e)
                    foo = 'bar'


            self.done_list.add(self.current_remote_dir)
            if self.links_to_visit_current:
                self.current_remote_dir = os.path.join(self.current_remote_dir, self.links_to_visit_current[0])
            else:
                if self.last_fork_location:
                    self.current_remote_dir = self.last_fork_location
                    self.last_fork_location = None
                else:
                    while True:
                        self.current_remote_dir = os.path.dirname(self.current_remote_dir)
                        if self.current_remote_dir == self.remote_base_dir or self.current_remote_dir + '/' == self.remote_base_dir:
                            if self.current_remote_dir not in self.done_and_empty_list and self.current_remote_dir + '/' not in self.done_and_empty_list:
                                self.walking_complete = False
                                break
                            else:
                                self.walking_complete = True
                                break
                        if self.current_remote_dir not in self.done_and_empty_list and self.current_remote_dir + '/' not in self.done_and_empty_list:
                            self.walking_complete = False
                            break
                    if self.walking_complete:
                        break
        
        return self.error_df_lists


    def _check_for_replicates(self):
        # Then we can count how many there are, add current dir to done
        # and continue walking the directories
        if len(self.fastq_gz_list_current) > 2:
            # For the sediments there are multiple sample IDs worth of fastq
            # files in a single directory so we need to check for the ratio
            # of fastq files to sample IDs.
            try:
                sample_id_set = set(['_'.join(_.split('_')[:2]) for _ in self.fastq_gz_list_current])
            except IndexError:
                foo = 'bar'
            if len(self.fastq_gz_list_current) > len(sample_id_set) * 2:
                # This is where we need to categorise the replication as
                # one of three different categories,
                # When all seq files have the same base_name, then we will call
                # this 'sequencing_replication'.
                # Where the base_name is different, but the pcr and dna names are the
                # same, we will call this 'unknown_replication'
                # Finally where there are any differences in pcr or dna names
                # I will call these 'method_replication'
                try:
                    base_names = {'-'.join(_.split('-')[:-1]) for _ in self.fastq_gz_list_current}
                    if '' in base_names:
                        print(self.fastq_gz_list_current)
                        foo = 'bar'

                except IndexError:
                    foo = 'bar'
                if len(base_names) > 1:
                    try:
                        self._process_unkn_method_replication(base_names, sample_id_set)
                    except IndexError:
                        foo = 'bar'
                else:
                    try:
                        self._process_seq_replication()
                    except IndexError:
                        foo = 'bar'

    def _process_unkn_method_replication(self, base_names, sample_id_set):
        # Then we need to check their PCR_sample_name and DNA_sample_name
        # To see if they are different.
        # This will produce two classes of seq difference
        # One where ReadSet names are different, but PCR and DNA are the same
        # Have to ask why?
        # The other set will be where ReadSet names are different, and one
        # of the PCR or DNA names are different.
        pcr_names_set = set()
        dna_names_set = set()
        temp_error_list = []
        for i, base_name in enumerate(base_names):
            try:
                sample_id = base_name.split('_')[1]
            except IndexError:
                foo = 'bar'
            try:
                if 'BID' in base_name:
                    element_one = base_name.split('-')[-2]
                    element_two = base_name.split('-')[-3].split('_')[-1]
                else:
                    element_one = base_name.split('-')[-1]
                    element_two = base_name.split('-')[-2].split('_')[-1]
            except IndexError:
                foo = 'bar'
            for j in range(10):
                readset_str = f'{j}_{element_two}.{element_one}'
                index_list = [_ for _ in self.readset_df[self.readset_df['sample_id'] == sample_id].index if readset_str in _]
                if len(index_list) == 1:
                    index_to_use = index_list[0]
                    break
                if len(index_list) > 1:
                    bid_name = base_name.split('-')[-1]
                    new_index_list = [_ for _ in index_list if bid_name in _]
                    if len(new_index_list) == 1:
                        index_to_use = new_index_list[0]
                        break
                if j == 9:
                    raise RuntimeError
            # if self.readset_df.at[index_to_use, 'sample_id'] != sample_id:
            #     raise RuntimeError('sample ids do not match')
            pcr_sample_name = self.readset_df.at[index_to_use, 'pcr_sample_name']
            dna_sample_name = self.readset_df.at[index_to_use, 'dna_sample_name']
            pcr_fl_sample_name = self.readset_df.at[index_to_use, 'pcr_fl_sample_name']

            pcr_names_set.add(pcr_sample_name)
            dna_names_set.add(dna_sample_name)
            fwd_read = [_ for _ in self.fastq_gz_list_current if (base_name in _) and ('R1' in _)][0]
            temp_error_list.append([sample_id, index_to_use, fwd_read, fwd_read.replace('R1', 'R2'), pcr_sample_name, pcr_fl_sample_name, dna_sample_name])
        # print(fastq_gz_list)
        if (len(pcr_names_set) != len(dna_names_set)) or (len(pcr_names_set) > len(sample_id_set)):
            self._log_method_replication(temp_error_list)
        else:
            self._log_unknown_replication(temp_error_list)

    def _log_unknown_replication(self, temp_error_list):
        # Then this a unknown_replication
        for temp_list in temp_error_list:
            another_temp_list = []
            for item in temp_list:
                another_temp_list.append(item)
            another_temp_list.extend(['unknown_replication', 'yellow', self.current_remote_dir])
            self.error_df_lists.append(another_temp_list)
    
    def _log_method_replication(self, temp_error_list):
        # then this is a 'method_replication'
        for temp_list in temp_error_list:
            another_temp_list = []
            for item in temp_list:
                another_temp_list.append(item)
            another_temp_list.extend(['method_replication', 'red', self.current_remote_dir])
            self.error_df_lists.append(another_temp_list)

    def _process_seq_replication(self):
        # Then this is a case of sequence_replication
        # We should be able to find a readset per fastq.
        # The readset should contain two bits of information
        # in the fastq and it should containing the -1 or -2
        # This is a pain in the arse!
        if len(self.fastq_gz_list_current) != 4:
            raise NotImplementedError

        # Do for read one
        read_1 = [_ for _ in self.fastq_gz_list_current if '-1_R1' in _][0]
        read_2 = [_ for _ in self.fastq_gz_list_current if '-2_R1' in _][0]
        for i, read_num in enumerate([read_1, read_2]):
            element_one = read_num.split('-')[-2]
            element_two = read_num.split('-')[-3].split('_')[-1]
            readset_str = f'{i + 1}_{element_two}.{element_one}'
            sample_id = read_num.split('_')[1]
            # now we need to look for this string in the indices
            # if we find more than one then rais error
            index_list = [_ for _ in self.readset_df.index if readset_str in _]
            if len(index_list) != 1:
                raise RuntimeError
            else:
                index_to_use = index_list[0]

            pcr_sample_name = self.readset_df.at[index_to_use, 'pcr_sample_name']
            dna_sample_name = self.readset_df.at[index_to_use, 'dna_sample_name']
            pcr_fl_sample_name = self.readset_df.at[index_to_use, 'pcr_fl_sample_name']
            self.error_df_lists.append([sample_id, index_to_use, read_num, read_num.replace('R1', 'R2'), pcr_sample_name, pcr_fl_sample_name, dna_sample_name, 'sequencing_replicate', 'green', self.current_remote_dir])

    def _make_auth_tup(self):
        auth_path = os.path.join(self.exe_path, 'auth.txt')
        with open(auth_path, 'r') as f:
            auth_lines = [line.rstrip() for line in f]
        return (auth_lines[0], auth_lines[1])
    
    def _make_readset_info_dir(self):
        # read in the three sepearate csv files
        coral_readset_df = pd.read_csv(os.path.join(self.readset_info_dir, "coral_readset_info.csv"), skiprows=[0], names=['readset', 'primers', 'sample_id', 'pcr_sample_name', 'dna_sample_name'])
        sed_readset_df = pd.read_csv(os.path.join(self.readset_info_dir, "ssed_readset_info.csv"), skiprows=[0], names=['readset', 'primers', 'sample_id', 'pcr_sample_name', 'dna_sample_name'])
        fish_readset_df = pd.read_csv(os.path.join(self.readset_info_dir, "fish_readset_info.csv"))
        # fish_readset_df.drop(columns='PCR FL sample name', inplace=True)
        fish_readset_df.columns = ['readset', 'primers', 'sample_id', 'pcr_sample_name', 'pcr_fl_sample_name', 'dna_sample_name']
        plankton_readset_df = pd.read_csv(os.path.join(self.readset_info_dir, "plankton_readset_info.csv"))
        # plankton_readset_df.drop(columns='PCR FL sample name', inplace=True)
        plankton_readset_df.columns = ['readset', 'primers', 'sample_id', 'pcr_sample_name', 'pcr_fl_sample_name', 'dna_sample_name']
        df = pd.concat([coral_readset_df, sed_readset_df, fish_readset_df, plankton_readset_df])
        return df.set_index('readset', drop=True)

# ReplicationWalkerHandler(marker='its2').start_rep_walker_mp()
# ReplicationWalkerHandler(marker='18s').start_rep_walker_mp()
ReplicationWalkerHandler(marker='16s_45').start_rep_walker_mp()
ReplicationWalkerHandler(marker='16s_full_45').start_rep_walker_mp()