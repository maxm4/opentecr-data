from pandas import read_csv, DataFrame, read_excel
from openpyxl import load_workbook
from itertools import islice
from numpy import nan, diff
from warnings import warn
import sigfig
import json, re, os

# homogenize the charge format    
def charge_format(master_reaction):       
    remove_string = None
    if re.search('\w(\d\-)', master_reaction):
        remove_string = re.search('\w(\d\-)', master_reaction).group(1)
        master_reaction = re.sub(remove_string, '-{}'.format(remove_string), master_reaction)
    return master_reaction, remove_string

# add the units of logarithm to the Magnesium concentration
def isnumber(string):
    if string != 'nan' and string is not nan:
        try:
            float(string)
        except:
            try:
                int(string)
            except:
                return False 
        return True

class merge_TECRs():
    def __init__(self, master_path, verbose = False):
        self.master_path = master_path
        self.verbose = verbose
        if not os.path.exists('TECRDB'):
            os.mkdir('TECRDB')
        with open(master_path, encoding='utf-8') as master:
            self.master_file = read_csv(master, dtype = 'object')
            self.master_file.fillna(' ')
            self.master_file.index.name = 'NIST_index'
            del self.master_file['Unnamed: 0']
            
            # insert the relevant columns
            empty_col = [' ' for row in range(len(self.master_file))]  
            self.master_file.insert(0, 'noor_index', empty_col)  
            self.master_file.insert(0, 'du_index', empty_col)
            self.master_file.insert(0, 'NIST_index', self.master_file.index)
            self.master_file.insert(4, 'KEGG Reaction:', empty_col)
            self.master_file.insert(5, 'CID Reaction:', empty_col)
        
    def merge(self, new_path, new_enzyme_col, new_reference_col, manual_curation_csv_path, excel_sheet):
        self.scraping = new_path.split('_')[-1].split('.')[0]
        print(self.scraping)
        if '.csv' in new_path:
            with open(new_path, encoding='utf-8') as new:
                self.new_file = read_csv(new, dtype = 'object')
                self.new_file.fillna(' ')
        elif '.xlsx' in new_path:
            ws = load_workbook(new_path)[excel_sheet]
            data = (islice(r, 1, None) for r in list(ws.values))
            self.new_file = DataFrame(data, index=[r[0] for r in ws.values], columns=next(ws.values)[1:])
            self.new_file.drop(self.new_file.index[0], inplace=True)
            self.new_file.insert(0, 'Enzyme', self.new_file.index)
            self.new_file.index = [x for x in range(len(self.new_file.index))]
            self.new_file.fillna(' ')
            # print(self.new_file.head())
        self._add_rows(new_enzyme_col, new_reference_col)
        self._merge_existing()
        self._incorporate_manual_curation(manual_curation_csv_path)
        self._confirm_merging()
        
        # export the processes database file
        self.master_file.to_csv(os.path.join('TECRDB', 'merged_TECRDB.csv'))
        # sleep(2)
        # with ZipFile('TECRDB.zip', 'w', compression = ZIP_LZMA) as zip:
        #     for file in ['TECR_consolidated.json', 'amalgamated_TECR_scrape.csv', 'TECRDB_scrape.csv']:
        #         zip.write(file)
        #         os.remove(file)
        return self.master_file
            
    def _merge_existing(self):
        def equal(val1, val2):
            if isinstance(val1, str) and isinstance(val2, str):
                if isnumber(val1) and isnumber(val2):
                    if self.__rounding(val1) != self.__rounding(val2):
                        error = r''.join([str(x) for x in [master_index, '___', val1, '___', val2]])
                        if self.verbose:
                           warn(error) 
                        errors.append(error)
                elif val1 != ' ' and val2 != ' ':
                    return val1 == val2
            else:
                return self.__rounding(val1) == self.__rounding(val2)
                
        matched_master_indices, errors_dictionary = {}, {}
        unmatched_entries = 0
        for new_index, new_row in self.new_file.iterrows():
            if new_index not in self.new_additions:
                matched_datum = False
                
                # define the new values for this datum 
                (NIST_index, du_index, noor_index, new_enzyme, add_kegg_reaction, add_cid_reaction, add_ec, new_reaction, reference_string, 
                 new_reference, new_temperature, new_ph, new_k, km, environmental_conditions, ionic_strength_1, enthalpy, solutes_1, buffer, add_pmg, 
                 add_method, add_ionic_strength, solutes_2 
                 ) = self.__define_row(new_index, new_row)

                # determine the set of possible matches
                errors = []
                matching_master_subset = self.master_file.loc[self.master_file['Enzyme:'] == new_enzyme]    #  & (master_file['Keq'] == new_k) & (master_file['T [K]'] == new_temperature) & (master_file['pH '] == new_ph)   #!!! other criteria besides the name should be considered for matching.
                for master_index, master_row in matching_master_subset.iterrows():  
                    # remove previously matched rows
                    if master_index in matched_master_indices:
                        error = f'IndexAlreadyMatched: The NIST index {master_index} is already matched to the new_index {matched_master_indices[master_index]}.'
                        if self.verbose:
                            warn(error)
                        errors.append(error)
                        continue

                    # match references
                    if not equal(new_reference, re.sub('_.+', '', master_row['Reference ID:'])):
                        error = f'ReferenceAlreadyAssigned: The NIST index {master_index} is already matched with the {re.sub("_.+", "", master_row["Reference ID:"])} reference, and thus cannot be assigned with the {new_reference}.'
                        if self.verbose:
                            warn(error)
                        errors.append(error)
                        continue

                    # compare values
                    if re.search('\w(\?\w+)', str(master_row['Keq'])):
                        master_row['Keq'] = re.sub('(\?\w+)', '', str(master_row['Keq']))
                    for new, master in zip(
                            [str(new_k).strip('~?'), str(new_temperature).replace('l', '1'), str(new_ph).strip('?~')], 
                            [str(master_row['Keq']).strip('~?'), str(master_row['T [K]']).replace('l', '1'), str(master_row['pH']).strip('?~')]):
                        if not equal(new, master):
                            error = f'DatumNotMatchErorr: The new value {new} | {new_index} does not match the existing value {master} | {master_index}.'
                            if self.verbose:
                               warn(error)
                            errors.append(error)
                            continue
                        
                    # match the reactions after multiple iterations 
                    master_reaction = master_row['Reaction:'] 
                    if self.scraping == 'noor':    
                        master_reaction = re.sub('\u00ce\u00b1|\u00ce\u00b2', '', master_reaction)
                        master_reaction = re.sub('\u00cf\u2030', '-w', master_reaction)
                    tries = 0
                    match = True
                    while new_reaction != master_reaction:
                        if tries == 0:
                            if re.search('= -\w', master_reaction):
                                remove_string = re.search('=(\s-)\w', master_reaction).group(1)
                                master_reaction = re.sub(remove_string, '-', master_reaction)
                        elif tries == 1:
                            if new_reaction != master_reaction:
                                if re.search(' -D-', master_reaction):
                                    remove_string = re.search('\s(-D-)', master_reaction).group(1)
                                    master_reaction = re.sub(remove_string, 'D-', master_reaction, 1)
                        elif tries == 2:
                            if re.search('\w\d\-', master_reaction):
                                loop = True
                                while loop:
                                    master_reaction, remove_string = charge_format(master_reaction)
                                    if remove_string is None:
                                        loop = False
                        elif tries == 3:
                            if re.search('\(\w\)\-', master_reaction):
                                master_reaction = re.sub('\(\w\)\-', '', master_reaction)
                        elif tries == 4:
                            if re.search('-lipoate', master_reaction):
                                master_reaction = re.sub('(-lipoate)', 'lipoate', master_reaction, 1)
                        else:
                            error = f'CodeError: The master reaction {master_reaction} | {master_index} does not match the new reaction {new_reaction} | {new_index}.'
                            if self.verbose:
                               warn(error)
                            errors.append(error)
                            match = False
                            break
                        tries += 1

                    # define the new data of the master file
                    if match:
                        matched_datum = True
                        matched_master_indices[master_index] = new_index
                        self.__redefine_master(master_row, new_row, master_index, new_index)
                        break

                if not matched_datum:
                    if self.verbose:
                        warn(f'CodeError: Failed index to match {new_enzyme} | {new_index}')
                    unmatched_entries += 1
                    errors_dictionary[new_index] = errors

        # test for standard_id uniqueness and unmatched values
        duplicate_entries = '\t'.join(
            [master_index for master_index in matched_master_indices if list(matched_master_indices.keys()).count(master_index) > 1])
        print('repeated entries:\n', duplicate_entries) 
        print('Unmatched indices: ', unmatched_entries)

        # export the unmatched datums
        with open(f'TECRDB/unmatched_{self.scraping}_TECRDB_datums.json', 'w') as output:
            json.dump(errors_dictionary, output, indent = 3)

    def _incorporate_manual_curation(self, manual_curation_path):
        # import and clean the manual curation file
        print(self.scraping)
        with open(manual_curation_path, 'r') as curation:
            manual_curation_csv = read_csv(curation, header = 0)
        headings = []
        for column in manual_curation_csv:
            headings.append(column.strip())
            for index, val in manual_curation_csv[column].iteritems():
                manual_curation_csv[column].iloc[index] = str(val).strip()
        manual_curation_csv.columns = headings

        # parse the manually curated content
        parsing_errors = []
        self.duplicates = set()
        for index, row in manual_curation_csv.iterrows():
            # characterize the curated datums            
            error = row['Error resolution']
            add = merge = False
            new_ids = [row['New index']]
            if row['Master file index'] == 'New' or re.search('sigfig', error):
                add = True
            elif row['Master file index'] != '--':   
                if all(isnumber(x) for x in row['Master file index'].split('-')):
                    master_file_ids = list(map(int, row['Master file index'].split('-')))
                else:
                    print(row['Master file index'].split('-'))
                merge = True
                if len(master_file_ids) > 1:
                    master_file_ids = [id_num + master_file_ids[0] for id_num in range(diff(master_file_ids)[0] + 1)]    
                    
            if all(isnumber(x) for x in new_ids[0].split('-')):
                new_ids = list(map(int, new_ids[0].split('-')))
            else:
                print(new_ids[0].split('-'))
            if len(new_ids) > 1:
                new_ids = [id + new_ids[0] for id in range(diff(new_ids)[0] + 1)]            
            
            # add new datums to the master file
            parsing_errors = None
            if add or merge:
                for new_id_index, new_id in enumerate(new_ids):
                    new_row = self.new_file.iloc[new_id]
                    if new_id not in self.master_file[self.scraping+'_index'].to_list():
                        if add:
                            self.master_file.loc[len(self.master_file)] = self.__define_row(new_id, new_row)
                        elif merge:
                            self.__redefine_master(self.master_file.iloc[new_id_index], new_row, master_file_ids[new_id_index], new_id)
                    else:
                        if self.verbose:
                            warn(f'CurationError: Repeated {self.scraping} index {new_id}.')
            else:
                # characterize the curated datums            
                if re.search('Duplicate', row['Error resolution']):
                    self.duplicates.update(new_ids)
                    print(f'The {new_ids} new_id is a duplicate.')
                elif re.search('already', row['Error resolution']):
                    print(f'The {new_ids} new_id is already matched.')
                else:
                    parsing_errors.append(master_file_ids)
                    print(f'CodeError: The {new_ids} new_id was not captured during parsing.')

        print('Parsing errors: ', parsing_errors)
    
    def _confirm_merging(self):
        missing_indices = set(self.new_file.index) - set(self.master_file[self.scraping+'_index'].to_list())
        missing_unique_indices = missing_indices - self.duplicates
        if missing_unique_indices != set():
            missing_unique_indices = sorted(missing_unique_indices)
            print(self.master_file[self.scraping+'_index'].to_list())
            print('missing unique indices', missing_unique_indices)
        else:
            print('all indices are captured')
            return missing_unique_indices
        
    def _add_rows(self, new_enzyme_col, new_reference_col, master_enzyme_column_name = 'Enzyme:', master_reference_col = 'Reference ID:'):
        # compare the enzymes
        new_enzymes = set(self.new_file[new_enzyme_col])  
        master_enzymes = set([enzyme.strip() for enzyme in self.master_file[master_enzyme_column_name].to_list() if enzyme != ' '])
        missing_master_enzymes = self.__set_contrast('enzymes', master_enzymes, new_enzymes, 'new file')
        
        # compare the references
        new_references = set(self.new_file[new_reference_col])
        master_references = set([re.sub('_.+', '', reference) for reference in self.master_file[master_reference_col].to_list() if reference != ' '])
        missing_master_references = self.__set_contrast('references', master_references, new_references, 'new file')
        original_master_length = len(self.master_file)
        
        # add new data rows
        self.new_additions = set()
        for new_index, new_row in self.new_file.iterrows():
            if new_row[new_reference_col] in missing_master_references or new_row[new_enzyme_col] in missing_master_enzymes:
                self.master_file.loc[len(self.master_file)] = self.__define_row(new_index, new_row)
                self.new_additions.add(new_index)
       
        # format the magnesium potential
        undescribed = list(self.master_file['Experimental conditions'])[original_master_length+1:]
        for index, row in enumerate(undescribed):
            if isnumber(row):
                self.master_file.at[index+original_master_length+1, 'Experimental conditions'] = f'{row} = -log[Mg+2]'

        if original_master_length == len(self.master_file):
            warn('CodeError: The master file has not changed length.')
        print('total additions', len(list(self.new_additions)))

    def __rounding(self, number):
        if self.scraping == 'noor':
            return float(number)
        elif self.scraping == 'du':
            return sigfig.round(number, 2)
        
    def __define_row(self, new_index, new_row):
        du_index = noor_index = kegg_reaction = cid_reaction = environmental_conditions = buffer = p_mg = None
        if self.scraping == 'noor':
            noor_index = new_index
            enzyme = new_row['enzyme_name']
            kegg_reaction = new_row['reaction']
            ec = new_row['EC']
            reaction = new_row['description']
            ref_id = new_row['reference']
            temperature = new_row['temperature']
            pH = new_row['p_h']
            keq = new_row['K_prime']
            if keq is nan:
                keq = new_row['K']
            ionic_strength = new_row['ionic_strength']
            p_mg = new_row['p_mg']
            method = new_row['method']
        elif self.scraping == 'du':
            environmental_conditions = ' ; '.join([str(x) for x in [new_row['media conditions'], new_row['electrolytes'], new_row['pMg']] if x is not None])
            environmental_conditions = re.sub(';\s+;', '', environmental_conditions)
            du_index = new_index
            enzyme = new_row['Enzyme']
            cid_reaction = new_row['Reaction formula in CID format']
            ec = new_row['EC value']
            reaction = new_row['Reaction']
            ref_id = new_row['Reference_id']
            temperature = new_row['T(K)']
            pH = new_row['pH']
            keq = new_row['K\'']
            ionic_strength = new_row['Ionic strength']
            buffer = new_row['Buffer/reagents/solute added']
            method = new_row['Method']
        return [
            None,                      # NIST_index 
            du_index,                  # f'https://w3id.org/related-to/doi.org/10.5281/zenodo.5494490/files/TableS1_Keq.csv#Keq_{new_index}',  
            noor_index,                # f'https://w3id.org/related-to/doi.org/10.5281/zenodo.3978439/files/TECRDB.csv#entry{new_index}', 
            enzyme,                    # add_enzyme
            kegg_reaction,             # add_kegg_reaction
            cid_reaction,              # add_cid_reaction
            ec,                        # add_ec 
            reaction,                  # add_reaction
            None,                      # reference_string
            ref_id,                    # add_reference
            temperature,               # add_temperature
            pH,                        # add_ph
            keq,                       # Keq  
            None,                      # km
            environmental_conditions,  # environmental_conditions
            ionic_strength,            # add_ionic_strength
            None,                      # enthalpy
            None,                      # solutes_1
            buffer,                    # buffer
            p_mg,                      # add_pmg
            method,                    # add_method
            None,                      # ionic_strength_1
            None,                      # solutes_2
            ]
    
    # merging values between matched datum
    def __redefine_master(self, master_row, new_row, master_index, new_index, verbose = False):
        # print the datum pair for manual inspection
        if verbose:
            print('\nmatched pair:', '\n', '='*len('matched pair:'))
            print('new_index', new_index)
            print('master_index', master_index, '\n')

        # match KEGG reactions
        master_kegg = master_row['KEGG Reaction:']
        if master_kegg == ' ' and self.scraping == 'noor':
            self.master_file.at[master_index, 'KEGG Reaction:'] = new_row['reaction']

        # match magnesium concentrations
        if self.scraping == 'noor':
            new_pmg = new_row['p_mg']
        elif self.scraping == 'du':
            new_pmg = new_row['pMg']
        if master_row['Experimental conditions'] != f'{new_pmg} = -log[Mg+2]' and new_pmg != ' ':
            if master_row['Experimental conditions'] in [None, ' ', nan]:
                self.master_file.at[master_index, 'Experimental conditions'] = f'{new_pmg} = -log[Mg+2]'
            else:
                self.master_file.at[master_index, 'Experimental conditions'] = master_row['Experimental conditions']+f" or {new_pmg} (-log[Mg+2])"
            if verbose:
                print(master_index, '\t', 'new pmg', '\t', self.master_file.at[master_index, 'Experimental conditions'])            
                
#         add_conditions = ' _ '.join([str(du_row['media conditions']), str(du_row['electrolytes']), ''])
#         add_conditions = re.sub('_\s+_', '', add_conditions)
#         if master_pmg in empty:
#             master_file.at[master_index, 'Experimental conditions'] = add_conditions
#             print(master_index, '\t', 'new pmg', '\t', master_file.at[master_index, 'Experimental conditions'])
#         else:
#             master_file.at[master_index, 'Experimental conditions'] = ' _ '.join([str(master_pmg), str(add_conditions)]) 
#             print(master_index, '\t', 'new pmg', '\t', master_file.at[master_index, 'Experimental conditions'])  
                
        # match methods
        if self.scraping == 'noor':
            new_method = new_row['method']
        elif self.scraping == 'du':
            new_method = new_row['Method']
        if master_row['Method:'] in [None, ' ', nan]:
            self.master_file.at[master_index, 'Method:'] = new_method
        elif new_method not in [None, ' ', nan]:
            new_methods = new_method.split(' and ')
            master_methods = master_row['Method:'].split(' and ')
            if set(new_methods + master_methods) != set():
                self.master_file.at[master_index, 'Method:'] = ' and '.join([method for method in set(new_methods + master_methods)])
            if verbose:
                print(f'The master index {master_index} possesses a new method: {self.master_file.at[master_index, "Method:"]}')

        # match EC values
        master_ec = master_row['EC Value:']
        if self.scraping == 'noor':
            new_ec = new_row['EC']
        elif self.scraping == 'du':
            new_ec = new_row['EC value']
        if master_ec not in [None, ' ', nan]:
            if master_ec == new_ec:
                self.master_file.at[master_index, 'EC Value:'] = master_ec   
            elif new_ec not in [None, ' ', nan]:
                if re.search('&', new_ec):
                    new_ecs = new_ec.split('&')
                    total_ecs = set(new_ecs + list(master_ec))
                else:
                    total_ecs = [master_ec, new_ec]
                self.master_file.at[master_index, 'EC Value:'] = " & ".join(total_ecs)
                if verbose:
                    print(f'The master index {master_index} possesses a new EC: {self.master_file.at[master_index, "EC Value:"]}')
        else:
            self.master_file.at[master_index, 'EC Value:'] = new_ec

        # match ionic strength concentrations
        master_ionic_strength = master_row['Ionic strength [mol/dm^3]']
        if self.scraping == 'noor':
            new_ionic_strength = new_row['ionic_strength']
        elif self.scraping == 'du':
            new_ionic_strength = new_row['Ionic strength']
        if master_ionic_strength in [None, ' ', nan]:
            self.master_file.at[master_index, 'Ionic strength [mol/dm^3]'] = new_ionic_strength
        elif new_ionic_strength not in [None, ' ', nan]:
            self.master_file.at[master_index, 'Ionic strength [mol/dm^3]'] = " & ".join([str(master_ionic_strength), str(new_ionic_strength)])
            if verbose:
                print(master_index, '\t', 'new ionic strength', '\t', self.master_file.at[master_index, 'Ionic strength [mol/dm^3]'])            

        # match the standard_id
        if master_row[self.scraping+'_index'] != ' ':
            if self.verbose:
                warn('The master_index {master_index} is predefined as {master_row[self.scraping+"_index"]}.')
            self.master_file.loc[len(self.master_file)] = self.__define_row(new_index, new_row)
        else:
            self.master_file.at[master_index, self.scraping+'_index'] = new_index
    
    # define the printing function
    def __set_contrast(self, data_description, master_set, set_2, set_2_description):
        # print the original sets
        print('='*30)
        print(f'\n{data_description} in the master file: ', len(master_set))
        print(f'{data_description} in the new file: ', len(set_2))

        # contrast the sets
        print(f'\nExtra {data_description} in the master file, versus {set_2_description}: ', len(master_set - set_2))
        print(f'\nMissing {data_description} in the master file, versus {set_2_description}: ', len(set_2 - master_set))
        if self.verbose:
            print('\n', master_set - set_2)
            print('\n', set_2 - master_set)
        return set_2 - master_set