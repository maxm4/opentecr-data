# import modules
from pandas import read_csv, isnull, DataFrame
from zipfile import ZipFile, ZIP_LZMA
from more_itertools import chunked
from pymed import PubMed
from io import StringIO
import httpx
import re, os

CHUNK_SIZE = 100

class mapRef():
    def __init__(self, mappings_url, zip_path = None, file_path = None, ):
        # load the CSV file and reference mappings
        request = httpx.request("GET", mappings_url)
        self.tecr_refs = read_csv(StringIO(request.content.decode("UTF-8")))
        if zip_path is not None:
            with ZipFile(zip_path, 'r') as zip:
                df = read_csv(zip.extract('amalgamated_TECR_scrape.csv'))
                os.remove('amalgamated_TECR_scrape.csv')
        elif file_path is not None:
            with open(file_path, 'r') as file:
                df = read_csv(file)     
        else:
            raise TypeError('The TECR data must be specified to complete mapping.')
        df = df.fillna(' ') # prevents spill-over of text
        df = df.astype(str)
        for column in df:
            if re.search('Unnamed', column):
                del df[column]
            # if column =='Reference ID:':
            #     df[column] = [re.sub('(_.+)', '', str(entry)) for entry in self.tecrdb[df]]
        self.tecrdb = df
        
    def apply(self, reference_ids_column, doi_column, pmid_column, export = False):
        # load references from PubMed
        tecr_refs_with_pubmed_id = self.tecr_refs[~isnull(self.tecr_refs.pmid)].copy()  # NOT null entries
        tecr_refs_with_pubmed_id["pmid"] = tecr_refs_with_pubmed_id.pmid.astype(int).astype(str) # remove decimals
        print(f"Collected {tecr_refs_with_pubmed_id.shape[0]} PubMed IDs")

        # parse the references for DOIs
        pubmed = PubMed(tool="MyTool", email="elad.noor@weizmann.ac.il")
        data = []
        count = 0
        for rows in chunked(tecr_refs_with_pubmed_id.itertuples(), 100):
            print(f'{count}/{tecr_refs_with_pubmed_id.shape[0]} references are mapped.', end='\r')
            pubmed_ids = " ".join([str(row.pmid) for row in rows])
            results = pubmed.query(pubmed_ids)
            for paper in results:
                try:
                    doi = paper.pubmed_doi
                except AttributeError:
                    doi = None

                pmid = str(paper.pubmed_id.split("\n")[0])
                authors = ", ".join([author["lastname"] + (" " + author["firstname"] if author["firstname"] else "") for author in paper.authors])
                data.append((pmid, doi, paper.publication_date.year, authors, paper.abstract))
            count += len(rows)

        # export the parsed information into a new CSV
        mappings = DataFrame(data=data, columns=["pmid", "doi", "year", "authors", "abstract"])
        self.mappings = mappings.join(tecr_refs_with_pubmed_id.set_index("pmid"), on="pmid", lsuffix="_from_pubmed", rsuffix="_from_Robert")
        
        # DOI and PMID columns are added to the TECRDB file
        new_column = [' ' for row in range(len(self.tecrdb))]
        self.tecrdb.insert(4, 'PMID', new_column)
        self.tecrdb.insert(5, 'DOI', new_column)

        references_added = 0
        for index, row in self.mappings.iterrows():
            if (row['pmid'] or row[doi_column]) != ' ':
                matching_master_subset = self.tecrdb.loc[self.tecrdb['Reference ID:'] == row[reference_ids_column]]
                for master_index, match in matching_master_subset.iterrows():
                    self.tecrdb.at[master_index, 'PMID'] = self.mappings.at[index, pmid_column]
                    self.tecrdb.at[master_index, 'DOI'] = self.mappings.at[index, doi_column]
                    references_added += 1
        print(f'References added to {references_added} datums')

        # export the dataframe
        self.tecrdb.to_csv('mapped_TECRDB_scrape.csv')
        self.mappings.to_csv("references_with_abstracts.csv")
        with ZipFile('TECRDB.zip', 'a', compression = ZIP_LZMA) as zip:
            for file in ['mapped_TECRDB_scrape.csv', "references_with_abstracts.csv"]:
                zip.write(file)
                os.remove(file)
            
        return self.tecrdb