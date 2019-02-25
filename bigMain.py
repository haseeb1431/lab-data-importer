import json

from datetime import datetime
from Bio import Entrez
from utils import config, bioweb_parser
from pg import pg_dbhelper
from time import sleep
from dateutil import parser

def search(query):
    Entrez.email = config.ncbi_db
    handle = Entrez.esearch(db=config.ncbi_db,
                            sort='relevance',
                            retmax='500',
                            retmode='xml',
                            term=query)
    results = Entrez.read(handle)
    return results


def search_author(query):
    Entrez.email = config.email
    handle = Entrez.esearch(db=config.ncbi_db,
                            sort='relevance',
                            retmax='1000',
                            retmode='xml',
                            term=query)
    results = Entrez.read(handle)
    return results


def fetch_details(id_list):
    ids = ','.join(id_list)
    Entrez.email = config.email
    handle = Entrez.efetch(db=config.ncbi_db,
                           retmode='xml',
                           id=ids)
    results = Entrez.read(handle)
    return results


def process_details():
    papers = fetch_details(['30715654'])
    for i, paper in enumerate(papers['PubmedArticle']):

        medlineCitation = paper['MedlineCitation']
        article = {
            'PMID': str(medlineCitation['PMID']),
            'ArticleTitle': medlineCitation['Article']['ArticleTitle'],
            'AbstractText': medlineCitation['Article']['Abstract']['AbstractText'][0],
            'JournalTitle': medlineCitation['Article']['Journal']['Title']
        }

        try:
            article['PublicationType'] = bioweb_parser. \
                parse_publication_type(medlineCitation['Article']['PublicationTypeList'])

            # article['Keywords'] = bioweb_parser.parse_keywords(medlineCitation['KeywordList'])
        except Exception as ex:
            print(ex)

        # article['Keywords'] = bioweb_parser.parse_keywords(medlineCitation['KeywordList'])
        print("%d) %s" % (i + 1, paper['MedlineCitation']['Article']['ArticleTitle']))


def fetch_summary(id_list):
    ids = ','.join(id_list)
    Entrez.email = config.email
    handle = Entrez.esummary(db=config.ncbi_db,
                             retmode='xml',
                             id=ids)
    results = Entrez.read(handle)
    return results


def process_summary(results):
    # TODO: Name of journal, list of pi, (last author), publication type
    results = results[0]

    # TODO handle special characters in name and title strings

    summary = {
        'FullJournalName': results['FullJournalName'],
        'PubTypeList': ",".join(results['PubTypeList']),
        'LastAuthor': results['LastAuthor'],
        'AuthorList': ",".join(results['AuthorList']),

        'Id': results["Id"],
        'Title': results["Title"]
    }

    if results["EPubDate"]:
        summary['PubDate'] = parser.parse(results["EPubDate"])
    elif results["PubDate"]:
        summary['PubDate'] = parser.parse(results["PubDate"])
    else:
        summary['PubDate'] = ""
    return summary


def getQueryTerm(author=None):
    # goes five year back from today
    query = '("{0}/{1}/{2}"[PDat]: "{3}/{1}/{2}"[PDat])'.format(datetime.today().year - 5,
                                                                datetime.today().month,
                                                                datetime.today().day,
                                                                datetime.today().year)
    if author is not None:
        query += '{0} AND {1}[Author]'.format(query, author)

    return query


if __name__ == '__main__':
    # Get all the papers from the author
    # we will use full name to avoid ambiguity among different authors
    authors = ['Stein Ulrike', 'Almut Nebel']

    for author in authors:
        query = getQueryTerm(author)
        search_results = search(query)
        # get the paper ids and loop over it
        id_list = search_results['IdList']
        print("found {0} Ids for query : {1}".format(len(id_list), query))
        # fetch details of all the ids
        for publication_id in id_list:
            try:
                summary_res = fetch_summary([publication_id])  # Optimize to send in batches
                paper_summary = process_summary(summary_res)
                paper_summary['PI'] = author

                # insert into the database
                insertSQL = pg_dbhelper.generateSQL([paper_summary])
                db_result = pg_dbhelper.insertData(insertSQL)
                print("{0} inserted into database".format(db_result))
            except Exception as ex:
                print(ex)

        # sleep it for a second to avoid the API limit
        sleep(1)

    # Pretty print the first paper in full to observe its structure
    print(json.dumps(paper_summary, indent=2, separators=(',', ':')))
