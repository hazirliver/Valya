from Bio.Entrez.Parser import StringElement
from Bio import Entrez
import pandas as pd
from typing import Dict, Union, List
from tqdm import tqdm
import time

Entrez.email = "Your.Name.Here@example.org"


def bacteria_species_to_id(bac_name: str) -> Dict[str, Union[str, None]]:
    """
    It takes a bacteria name as a string, searches the NCBI taxonomy database for the bacteria name, and returns the
    bacteria's taxonomy ID as a string

    :param bac_name: the name of the bacteria species you want to find the ID for
    :type bac_name: str
    :return: A dictionary with the bacteria name as the key and the taxonomy ID as the value.
    """
    handle = Entrez.esearch(db="taxonomy", retmax=1, term=bac_name)
    records = Entrez.read(handle)
    handle.close()
    return {bac_name: records['IdList'][0]}


def fetch_bacteria_lineage(bac_id: str) -> StringElement:
    """
    It takes a bacteria ID and returns the lineage of that bacteria

    :param bac_id: The ID of the bacteria you want to get the lineage of
    :type bac_id: str
    :return: A list of dictionaries.
    """
    handle = Entrez.efetch(db="taxonomy", id=bac_id)
    record = Entrez.read(handle)
    handle.close()

    return record[0]["LineageEx"]


def id_list_to_lineage(id_list: List[Dict]) -> Dict:
    """
    > This function takes a list of dictionaries, where each dictionary contains a bacteria name and a bacteria ID, and
    returns a dictionary where each bacteria name is mapped to a list of its lineage

    :param id_list: a list of dictionaries, where each dictionary has a bacteria name as the key and a bacteria ID as the
    value
    :type id_list: StringElement
    :return: A dictionary with the bacteria name as the key and the lineage as the value.
    """
    lineage_dict = {}
    for bacteria_dict in tqdm(id_list):
        for bac_name, bac_id in bacteria_dict.items():
            lineage_dict[bac_name] = fetch_bacteria_lineage(bac_id)
            time.sleep(1)
    return lineage_dict


def lineage_to_df(lineage_dict: Dict) -> pd.DataFrame:
    """
    It takes a dictionary of dictionaries and turns it into a dataframe

    :param lineage_dict: A dictionary of lineages, where the keys are the species names and the values are the lineages
    :type lineage_dict: Dict
    :return: A dataframe with the bacteria species as the index and the columns are the taxonomic levels.
    """
    df_dict = {}
    for bacteria_species, bacteria_lineage in lineage_dict.items():
        tmp = pd.DataFrame(lineage_dict[bacteria_species])
        tmp['species'] = [bacteria_species] * tmp.shape[0]
        df_dict[bacteria_species] = tmp

    return pd.concat(df_dict, axis=1).stack().T


def main() -> None:
    df = pd.read_csv('../Input/Bacteria_encoding.csv')
    df.replace({'Pseudomona aeruginosa': 'Pseudomonas aeruginosa'}, inplace=True)

    id_list = [bacteria_species_to_id(bacteria_species) for bacteria_species in tqdm(df.iloc[:, 0].tolist())]
    lineage_dict = id_list_to_lineage(id_list)
    tax_df = lineage_to_df(lineage_dict)
    tax_df.to_csv('../Output/tax_df.csv')


if __name__ == '__main__':
    main()
