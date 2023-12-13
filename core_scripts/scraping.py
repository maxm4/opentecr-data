# -*- coding: utf-8 -*-
"""
Created on Mon May 16 09:24:58 2022

@author: Andrew Freiburger
"""
from requests.packages.urllib3.util.retry import Retry
from requests.adapters import HTTPAdapter
from pandas import DataFrame, read_csv, concat
from zipfile import ZipFile #, ZIP_LZMA
from bs4 import BeautifulSoup
from warnings import warn
from shutil import move
from time import sleep
import requests
import numpy
import math, json, os, re

class TECRDB():
    def __init__(self, printing = True):
        # defining the website
        self.printing = printing
        
    def scrape(self,):
        root_url = "https://randr.nist.gov/enzyme/DataDetails.aspx?ID="
        end_url = "&finalterm=&data=enzyme"
        
        # identify the table and rows of pertinent data
        with open(os.path.join(os.path.dirname(__file__), 'Enzyme Thermodynamic Database.html'), 'r') as tecr_home_page:
            bs = BeautifulSoup(tecr_home_page, 'lxml')
        body = bs.find("table", attrs = {'id': 'MainBody_gvSearch'}).find_all("tr")
        total_entries = math.floor(1*len(body))

        # defining the boundaries of the dataframe section
        # open the requests session
        session = requests.Session()
        adapter = HTTPAdapter(max_retries = Retry(connect=3, backoff_factor=0.5))
        session.mount('https://', adapter)

        # loop through the enzyme id values 
        if not os.path.exists('TECR_scraping'):
            os.mkdir('TECR_scraping')

        index_count = loop_count = name_iteration = 0
        entry_dfs = []
        for id_row in range(1, total_entries):   
            id_value = body[id_row].find("a").text
            total_url = root_url + id_value + end_url
            bs = BeautifulSoup(session.get(total_url).text, 'lxml')
            
            # parsing the reaction names and strings   
            enzyme_name = body[id_row].find('span', attrs = {'id': 'MainBody_gvSearch_lblEnzyme_%s' %(name_iteration)}).text
            reaction = body[id_row].find('span', attrs = {'id': 'MainBody_gvSearch_lblReaction_%s' %(name_iteration)}).text
            name_iteration += 1
            
            ## first set of information 
            tables = bs.find_all("table", attrs={"id": "MainBody_DataList1"})
            if len(tables) == 0:
                warn(f'TECRError: The {id_value} reference {total_url} does not possess data.')
                continue
            body2 = tables[0].find_all("tr")
            body_rows2 = body2[1:]        
            information_entries_list, information_values_list, each_row2 = [], [], []
            for row in range(len(body_rows2)):
                for row_element in body_rows2[row].find_all("td"):
                    row_refined2 = re.sub("(\xa0)|(\n)|,","",row_element.text)
                    each_row2.append(row_refined2)

            column_count = 0
            for i, element in enumerate(each_row2):
                if i == 0 or i % 2 == 0:
                    information_entries_list.append(element)
                    column_count += 1
                else:
                    information_values_list.append(element)
                    column_count += 1
            column_count /= 2  #!!! why is the column count halved?
            
            table_df2 = DataFrame(
                data = [information_values_list], columns = information_entries_list, index = range(index_count, index_count+1)
                )
            table_df2.drop([table_df2.columns[-2], table_df2.columns[-1]], axis=1, inplace=True)
            
            ## second set of information 
            tables1 = bs.find_all("table", attrs = {"id": "MainBody_extraData"})
            if len(tables1) != 1:
                warn(f'The {id_value} reference {total_url} possesses an unexpected data structure.')
                continue
            if self.printing:
                print(id_value, f'\t\t{id_row}/{total_entries} enzymes', 
                      f'\t\t{index_count} datums', '\t\t\t\t', end = '\r')
            body1 = tables1[0].find_all("tr")
            heads, body_rows1 = body1[0], body1[1:]
            headings = ['Enzyme:', 'Reaction:']
            for head in heads.find_all("th"):
                head = (head.text).rstrip("\n")
                headings.append(head)
        
            total_rows = []
            for row_number in range(len(body_rows1)):
                each_row = [enzyme_name, reaction]
                for row_element in body_rows1[row_number].find_all("td"):
                    row_refined = re.sub("(\xa0)|(\n)|,","",row_element.text)
                    each_row.append(row_refined)
                total_rows.append(each_row)
            
            table_df1 = DataFrame(
                data = total_rows, columns = headings, index = range(index_count, len(body_rows1)+index_count)
                )
            table_df1.drop(table_df1.columns[len(table_df1.columns)-1], axis=1, inplace=True)
            
            # merge current and existing dataframes 
            this_dataframe = table_df1.join(table_df2)
            this_dataframe.index.name = 'index'
            if loop_count == 0:
                old_dataframe = this_dataframe
                old_dataframe.index.name = 'index'
            elif loop_count > 0:
                common_columns = list(set(this_dataframe.columns) & set(old_dataframe.columns))  # intersection operator
                current_dataframe = old_dataframe.merge(this_dataframe, on = common_columns, how = 'outer')
                old_dataframe = current_dataframe  
                
            # export the dataframe 
            index_count += len(body_rows1)
            time_delay = max_referenes_per_csv = 0
            sleep(time_delay)
            if loop_count == max_referenes_per_csv:
                # id_value = re.sub('(/)', '-', id_value)
                # while os.path.exists(os.path.join('TECR_scraping', f'{id_value}_{output_loop}.csv')):
                #     output_loop += 1
                entry_dfs.append(old_dataframe)
                loop_count = 0 
            else:        
                loop_count += 1 
                
        # combine all of the dataframes
        combined_df = concat(entry_dfs)
        if self.printing:
            display(combined_df)
        
        # refine the dataframe 
        combined_df = combined_df.fillna(' ') # prevents spill-over of text
        middle_dataframe_columns = ['T(K)', 'pH ', 'K<sub>c</sub>\' ', 'δ<sub>r</sub>H\'<sup>o</sup>(kJ.mol<sup>-1</sup>)', 'Km\'']
        left_dataframe_columns = ['index', 'Enzyme:', 'EC Value:', 'Reaction:', 'Reference:', 'Reference ID:'] 
        right_dataframe_columns = list(set(combined_df.columns) - set(left_dataframe_columns) - set(middle_dataframe_columns))
        self.scraped_df = combined_df.reindex(
            columns = left_dataframe_columns + middle_dataframe_columns + right_dataframe_columns
            )
        
        # export the dataframe
        if not os.path.exists('TECRDB'):
            os.mkdir('TECRDB')
        self.scraped_df.to_csv('TECRDB/TECRDB_scrape.csv')
        # with ZipFile('TECRDB.zip', 'w', compression = ZIP_LZMA) as zip:
        #     zip.write('TECRDB_scrape.csv')
        #     os.remove('TECRDB_scrape.csv')

    def amalgamate(self, zip_path = None):
        def merge_cells(re_search, col_name, printed):
            if re.search(re_search, this_column):
                solute = re.search(re_search, this_column).group(1)
                if row[this_column] not in [' ', '?']:
                    if row[col_name] == ' ':
                        df.at[index, this_column] = str(row[this_column]) + ' ' + solute
                    if row[this_column] != ' ':
                        if row[col_name] != (row[this_column] or str(row[this_column]) + ' ' + solute):
                            df.at[index, this_column] = str(row[col_name])+' & '+str(row[this_column])+' '+solute
                if this_column !=  col_name:
                    if not printed:
                        print('combined\t', this_column)
                        printed = True
                    combined_columns.add(this_column)
                    
            return printed
        
        # import the scraped data
        if not os.path.exists('TECRDB'):
            os.mkdir('TECRDB')
        if zip_path == None:
            df = read_csv(os.path.join('TECRDB/TECRDB_scrape.csv'))
        elif zip_path:  
            with ZipFile(zip_path, 'r') as zip:
                df = read_csv(zip.extract('TECRDB_scrape.csv'))
        df = df.fillna(' ') # prevents spill-over of text
        df = df.astype(str)
        if os.path.exists('TECRDB_scrape.csv'):
            move('TECRDB_scrape.csv', 'TECRDB/TECRDB_scrape.csv')
        
        combined_columns = set()
        re_searches = {
            '(?<=c\()(\w+\d?\+?)(?<!,)': 'c(glycerol,mol dm<sup>-3</sup>)',
            '(?<=m\()(\w+\d?\+?)(?<!,)': 'm(MgCl2,mol.kg<sup>-1</sup>)',
            '(ë«|Î\´|δ)': 'δ<sub>r</sub>H(cal)/kJ mol<sup>-1</sup>)',
            '(I<sub>c)':'I<sub>c</sub>(mol dm<sup>-3</sup>)',
            '(?i)(?!Km\' )(^K)': 'K<sub>c</sub>\' ', 
            '(Km\')(?! )':'Km\' ',
            }
        print('\nColumns:\n', '='*len('Columns:'))
        for this_column in df:
            printed = False
            if 'index' in this_column:
                continue
            for index, row in df.iterrows():
                # copy entries across all rows
                if row['Reference ID:'] != ' ':
                    reference_row = row
                else:
                    for row_col, val in row.iteritems():
                        if val == ' ':
                            df.at[index, row_col] = reference_row[row_col]
                
                # combine similar columns
                for re_search, col_name in re_searches.items():
                    printed = merge_cells(re_search, col_name, printed)
                if this_column == 'EC Value:':
                    ecs = []
                    for ec in row[this_column].split('&'):
                        if not re.search('((\d+\.)+(\d+|\-)?)', row[this_column]):
                            df.at[index, this_column] = ' '
                            break
                        ecs.append(re.search('((\d+\.)+(\d+|\-)?)', row[this_column]).group().strip())
                    if df.at[index, this_column] != ' ':
                        df.at[index, this_column] = ' & '.join(ecs)
                elif this_column in ['buffer(mol dm<sup>-3</sup>)', 'buffer and/or salt ', 'media ', 'buffer ']:
                    if not printed:
                        print('combined\t', this_column)
                        printed = True
                    if row[this_column] not in [' ', '?']:
                        if row['Buffer:'] == ' ':
                            df.at[index, 'Buffer:'] = row[this_column]
                        if row['Buffer:'] != ' ':
                            if not re.search(row[this_column], row['Buffer:']):
                                df.at[index, 'Buffer:'] = str(row['Buffer:']) + ' + ' + row[this_column]
                    combined_columns.add(this_column)
                elif this_column in ['salt ', 'cosolvent ', 'added solute ', 'protein ', 'added solute ', 
                                   'percent(dimethyl sulfoxide) ', 'p(MPa)']:
                    if not printed:
                        print('combined\t', this_column)
                        printed = True
                    if row[this_column] not in [' ', '?']:
                        if ['solvent '] == ' ':
                            df.at[index, 'solvent '] = row[this_column]
                            if this_column == 'p(MPa)':
                                df.at[index, 'solvent '] += ' MPa'  
                            elif this_column == 'percent(dimethyl sulfoxide) ':
                                df.at[index, 'solvent '] += ' % DMSO'   
                        else:
                            if not re.search(row[this_column], row['solvent ']):
                                df.at[index, 'solvent '] = row['solvent '] + '  +  ' + row[this_column]
                    combined_columns.add(this_column)
                else:
                    if not printed:
                        print('not combined\t', this_column)
                        printed = True
                        
        # delete the combined columns and export
        print('\nCombined columns:')
        for column in combined_columns:
            del df[column]
            print(column)
        for column in df.columns:
            if re.search('index|Unnamed', column):
                print(column)
                del df[column]
        df.rename(columns = {
            'I<sub>c</sub>(mol dm<sup>-3</sup>)':'Ionic strength [mol/dm^3]', 
            'δ<sub>r</sub>H(cal)/kJ mol<sup>-1</sup>)':'Enthalpy [kJ/mol]',
            'I<sub>m</sub>(mol.kg<sup>-1</sup>)':'Ionic strength [molal]', 
            'c(glycerol,mol dm<sup>-3</sup>)':'solutes [mol/dm^3]', 
            'm(MgCl2,mol.kg<sup>-1</sup>)':'solutes [molal]', 
            'solvent ':'Experimental conditions', 
            'K<sub>c</sub>\' ':'Keq', 
            'pMg ': '-log[Mg+2]',
            'Km\' ': 'Km\'',
            'T(K)':'T [K]',
            'pH ':'pH'
            }, inplace = True)
        
        self.amalgamated_df = df
        self.amalgamated_df.to_csv("TECRDB/amalgamated_TECRDB_scrape.csv") 
            
        # count down for processing and organizing the data
        def assign_values(col, reference_list, values_list, temperatures_list, ph_list, added):
            reference_list.append(False)   
            if row[col] != ' ':
                cleaned_keq = re.search('(\-?\d+\.?\d*)', str(row[col])).group()
                values_list.append(float(cleaned_keq)) 
                reference_list[-1] = True
                if not added:
                    temperatures.append(row['T [K]'])
                    ph_list.append('nan')
                    if row['pH'] != ' ':
                        ph_list[-1] = row['pH']
                    added = True
                        
            return reference_list, values_list, temperatures_list, ph_list, added
        
        # acquire a list of all enzymes
        enzyme_list = []
        for index, row in df.iterrows():
            if row['Enzyme:'] not in enzyme_list and row['Enzyme:'] != ' ':
                enzyme_list.append(row['Enzyme:'])            
        enzymes = [re.search('(\w.*)',original_enzyme).group() for original_enzyme in enzyme_list]
        data_per_enzyme = {}
        count = 1
        for enzyme in enzymes:
            print(f'The data is being assembled and organized ... {count}/425', end = '\r')
            
            # lists of the database variables
            added = False
            Keq_values, km_values, enthalpy_values, temperatures, phs= [], [], [], [], []
            references, reactions, kms, enthalpies, keqs = [], [], [], [], []
            for index, row in df.iterrows():
                iteration = 1
                if row['Enzyme:'] == ' '+enzyme:
                    reactions.append(row['Reaction:'])
                    references.append('Ibid')
                    if row['Reference:'] != ' ':
                        references[-1] = row['Reference:']
                    
                    # values are organized into lists
                    keqs, Keq_values, temperatures, phs, added = assign_values('Keq', keqs, Keq_values, temperatures, phs, added)                        
                    kms, km_values, temperatures, phs, added = assign_values('Km\'', kms, km_values, temperatures, phs, added)                        
                    enthalpies, enthalpy_values, temperatures, phs, added = assign_values(
                        'Enthalpy [kJ/mol]', enthalpies, enthalpy_values, temperatures, phs, added
                        )                        
                        
                    #loop through the unlabeled rows of each enzyme
                    while iteration + index < len(df) and df.at[index + iteration, 'Enzyme:'] == ' ':
                        keqs, Keq_values, temperatures, phs, added = assign_values('Keq', keqs, Keq_values, temperatures, phs, added)                        
                        kms, km_values, temperatures, phs, added = assign_values('Km\'', kms, km_values, temperatures, phs, added)                        
                        enthalpies, enthalpy_values, temperatures, phs, added = assign_values(
                            'Enthalpy [kJ/mol]', enthalpies, enthalpy_values, temperatures, phs, added
                            )         
                        iteration += 1
        
            # statistical processing of scraped values
            average_keq = std_keq = average_km = std_km = average_enthalpy = std_enthalpy = 'nan' 
            if len(Keq_values) != 0:
                average_keq, std_keq = numpy.mean(Keq_values), numpy.std(Keq_values)
            if len(km_values) != 0:
                average_km, std_km = numpy.mean(km_values), numpy.std(km_values)
            if len(enthalpy_values) != 0:
                average_enthalpy, std_enthalpy = numpy.mean(enthalpy_values), numpy.std(enthalpy_values) 
                
            #store the information into a nested dictionary structure
            data_per_enzyme[enzyme] = {'reaction':reactions,
                                       'experimental temperatures':temperatures,
                                       'experimental phs':phs,
                                       'keq reference':references,
                                       'Keq':{'keq values in the reference':keqs,
                                              'keqs':Keq_values, 
                                              'keq quantity':len(Keq_values), 
                                              'keq average':average_keq, 
                                              'keq standard deviation':std_keq
                                              },
                                       'Km':{'km values in the reference':kms,
                                             'km values':km_values,
                                             'km average':average_km,
                                             'km standard deviation':std_km
                                             },
                                       'Enthalpy':{'enthalpy values in the reference':enthalpies,
                                                   'enthalpy values':enthalpy_values,
                                                   'enthalpy average':average_enthalpy,
                                                   'enthalpy standard deviation':std_enthalpy
                                                   }
                                      }
            count += 1
        
        #export the database dictionary as a JSON file
        with open('TECRDB/TECRDB_consolidated.json', 'w') as output:
            json.dump(data_per_enzyme, output, indent = 4)
        # sleep(2)
        # with ZipFile('TECRDB.zip', 'w', compression = ZIP_LZMA) as zip:
        #     for file in ['TECR_consolidated.json', 'amalgamated_TECR_scrape.csv', 'TECRDB_scrape.csv']:
        #         zip.write(file)
        #         os.remove(file)
