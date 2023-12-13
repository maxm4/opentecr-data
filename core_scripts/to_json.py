# -*- coding: utf-8 -*-
"""
Created on Mon May 17 09:24:58 2022

@author: Andrew Freiburger
"""
from zipfile import ZipFile, ZIP_LZMA
from pandas import read_csv
from glob import glob
import json, os, re

class JSON_datum():
    def __init__(self, csv_path, openTECR_path):
        self.openTECR_path = openTECR_path
        self.csv = read_csv(csv_path).astype(str)
        
    def tecr_to_json(self, json_template_path, curator = 'APF'):
        with open(json_template_path, 'r') as template:
            json_template = json.load(template)
        for index, row in self.csv.iterrows():
            self.template = json_template
            for key, value in self.template.items():
                if key == 'curator':
                    if curator == 'APF':
                        value['orcid'] = 'https://orcid.org/0000-0002-7288-535X'
                        value['name'] = 'Andrew Philip Freiburger' 
                elif key == 'CuratedMeasurement':
                    # define datum content
                    value['reference']['pmid'] = row['PMID']
                    value['reference']['doi'] = row['DOI']
                    value['reference']['tecrdb_string'] = row['Reference:']
                    value['reference']['tecrdb_id'] = row['Reference ID:']
                    value['representative_reaction']['name'] = row['Enzyme:']
                    value['representative_reaction']['KEGG'] = row['KEGG Reaction:']
                    value['representative_reaction']['CID'] = row['CID Reaction:']

                    ## population the reaction content
                    if '=' not in row['Reaction:']:
                        value['representative_reaction']['stoichiometry'] = row['Reaction:']
                    else:
                        reactants, products = row['Reaction:'].split('=')
                        reactants_dict, products_dict = {}, {}
                        for reactant in reactants.split(' + '):
                            stoich = 1
                            if re.search('\d\s', reactant):
                                stoich = re.search('\d\s', reactant).group()
                            reactants_dict[reactant] = stoich
                        for product in products.split(' + '):
                            stoich = 1
                            if re.search('\d\s', product):
                                stoich = re.search('\d\s', product).group()
                            products_dict[product] = stoich

                        if value['representative_reaction']['stoichiometry'] != '':
                            value['representative_reaction']['stoichiometry']['reactants'] = reactants_dict
                            value['representative_reaction']['stoichiometry']['products'] = products_dict

                # define the measurement data
                if key == 'equilibriumConstant':
                    value = row['Keq']
                elif key == 'hydrogenPotential':
                    value = row['pH ']
                elif key == 'temperature':
                    value = row['T [K]']
                elif key == 'ionicStrength':
                    value = row['Ionic strength [mol/kg]']
                elif key == 'magnesiumPotential' and re.search('(\d+.?\d+) = -log\[Mg\+2\]', str(row['Experimental conditions'])):
                    value = re.search('(\d+.?\d+) = -log\[Mg\+2\]', str(row['Experimental conditions'])).group()

            self._export(row['Reference ID:'], '', 'tecr')
        self._export('NIST_TECRDB_JSONs', '', 'tecr')

    def half_rxn_to_json(self, json_template_path, curator = 'APF'):
        with open(json_template_path, 'r') as template:
            json_template = json.load(template)
        for index, row in self.csv.iterrows():
            self.template = json_template
            for key, value in self.template.items():
                if key == 'curation':
                    if curator == 'APF':
                        value['curator']['orcid'] = 'https://orcid.org/0000-0002-7288-535X'
                        value['curator']['name'] = 'Andrew Philip Freiburger'
                    value['comments'] = row['comment']
                elif key == 'CuratedMeasurement': 
                    # define datum content
                    value['reference']['doi'] = row['doi']
                    value['reference']['year'] = re.search('(?<=\()(\d+)(?=\))', row['reference']).group()
                    value['reference']['string'] = row['reference']
                    value['compound']['name'] = row['name']
                    value['hydrogenPotential'] = row['p_h']
                    value['temperature'] = row['temperature']
                    value['ionicStrength'] = row['ionic_strength']
                    value['standard_voltage_prime'] = row['standard_E_prime']
    #               value['equilibriumConstant'] = value['standard_voltage_prime']  #!!! The equilibrium constant must be calculated from the provided standard voltage
                    if 'metanetx' in row['CID_red']:
                        value['compound']['reduced']['MetaNetX']['id'] = row['CID_red']
                        value['compound']['reduced']['MetaNetX']['charge'] = row['charge_red']
                        value['compound']['reduced']['MetaNetX']['nH'] = row['nH_red']
                    if 'pubchem' in row['CID_red']:
                        value['compound']['reduced']['PubChem']['id'] = row['CID_red']
                        value['compound']['reduced']['PubChem']['charge'] = row['charge_red']
                        value['compound']['reduced']['PubChem']['nH'] = row['nH_red']
                    if 'metanetx' in row['CID_ox']:
                        value['compound']['oxidized']['MetaNetX']['id'] = row['CID_ox']
                        value['compound']['oxidized']['MetaNetX']['charge'] = row['charge_ox']
                        value['compound']['oxidized']['MetaNetX']['nH'] = row['nH_ox']
                    if 'pubchem' in row['CID_ox']:
                        value['compound']['oxidized']['PubChem']['id'] = row['CID_ox']
                        value['compound']['oxidized']['PubChem']['charge'] = row['charge_ox']
                        value['compound']['oxidized']['PubChem']['nH'] = row['nH_ox']

            self._export(row['name'], row['reference'], 'half_rxns')
        self._export('Elad_half_rxns', '', 'half_rxns')
                    
    def _export(self, name, reference, data_source, zip_contents = False):
        export_name = '_'.join([name, reference])
        if not zip_contents:
            # export the datum JSON
            if not os.path.exists(os.path.join(self.openTECR_path, 'datum_points', data_source)):
                os.mkdir(os.path.join(self.openTECR_path, 'datum_points', data_source))
            count = 0
            while os.path.exists(os.path.join(self.openTECR_path, 'datum_points', data_source, export_name+f'{count}.json')):
                count += 1
            with open(os.path.join(self.openTECR_path, 'datum_points', data_source, export_name+f'{count}.json'), 'w') as out:
                json.dump(self.template, out, indent = 4)
        else:
            # zip the set of JSON data
            with ZipFile(export_name+'.zip', 'a', compression = ZIP_LZMA) as _zip:
                for file in glob(os.path.join(self.openTECR_path, 'datum_points', data_source, '*.json')):
                    _zip.write(file)
                    os.remove(file)