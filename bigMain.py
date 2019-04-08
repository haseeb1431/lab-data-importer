import json

from datetime import datetime
from Bio import Entrez
from utils import config, bioweb_parser
from pg import pg_dbhelper
from time import sleep
from dateutil import parser

import azureDataPipeline

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
    """Generate query string which we will send to the API"""
    # goes five year back from today
    query = '("{0}/{1}/{2}"[PDat]: "{3}/{1}/{2}"[PDat])'.format(datetime.today().year - 5,
                                                                datetime.today().month,
                                                                datetime.today().day,
                                                                datetime.today().year)
    if author is not None:
        query = '{0} AND {1}[Author]'.format(query, author)

    return query


if __name__ == '__main__':
    # Get all the papers from the author
    # we will use full name to avoid ambiguity among different authors
    # 'Stein Ulrike'
    #authors = ['Almut Nebel']
    
    #return test_method();

    azureDataPipeline.run_pipeline()

    authors = [
         'Stephan Günther'
        , 'César Muñoz-Fontela'
        , 'Jonas Schmidt-Chanasit'
        , 'Maura Dandri'
        , 'Joachim Hauber'
        , 'Kay Grünewald'
        , 'Thomas Dobner'
        , 'Marcus Altfeld'
        , 'Wolfram Brune'
        , 'Gülsah Gabriel'
        , 'Adam Grundhoff'
        , 'Linlin Zhang'
        , 'Norbert Tautz'
        , 'Stefan Taube'
        , 'Susanne Herold'
        , 'Nina Babel'
        , 'Brigitte Biesinger-Zwosta'
        , 'Armin Ensser'
        , 'Thomas Gramberg'
        , 'Andrea Thoma-Kreß'
        , 'Michael Mach'
        , 'Manfred Marschall'
        , 'Frank Neipel'
        , 'Thomas Stamminger'
        , 'Linda Grosche'
        , 'Barbara Walch-Rückheim'
        , 'Susanne M. Bailer'
        , 'Holger Jeske'
        , 'Ralf Amann'
        , 'Barbara Adler'
        , 'Karl-Klaus Conzelmann'
        , 'Ulrike Protzer'
        , 'Volker Bruss'
        , 'Andreas Pichlmair'
        , 'Sabrina Schreiner'
        , 'Martina Anton'
        , 'Jennifer Altomonte'
        , 'Erdem Gültekin Tamgüney'
        , 'Ina Maja Vorberg'
        , 'Robert Daniels'
        , 'Anna-Lena Spetz'
        , 'Sonia Andersson'
        , 'Mariet C.W. Feltkamp'
        , 'Ann C.T.M. Vossen'
        , 'Martijn J. van Hemert'
        , 'Montserrat Bárcena'
        , 'Martijn J. van Hemert'
        , 'Manuel Goncalves'
        , 'Rob Hoeben'
        , 'Ronald van Rij'
        , 'Willem Melchers'
        , 'Anne Wensing'
        , 'Andy Hoepelman'
        , 'Emmanuel Wiertz'
        , 'Neeltje Kootstra'
        , 'Monique van Oers'
        , 'André Boonstra'
        , 'Qiuwei Pan'
        , 'Andrea Woltman'
        , 'Eric C.M. van Gorp'
        , 'Barry Rockx'
        , 'Georges M.G.M. Verjans'
        , 'Rob Gruters'
        , 'Ron Fouchier'
        , 'B.G. van den Hoogen'
        , 'Rik L. De Swart'
        , 'MPG Koopmans'
        , 'Mari Kaarbø'
        , 'Gunn Berit Olsson'
        , 'Katharina Lahl'
        , 'Graham Belsham'
        , 'Thomas Bruun Rasmussen'
        , 'Lars Erik Larsen'
        , 'Allan Randrup Thomsen'
        , 'Jens Bukh'
        , 'Dorte Kinggaard Holm'
        , 'Rasmus O. Bak'
        , 'Steffen Junker'
        , 'Yonglun Luo'
        , 'Bjørn Andersen Nexø'
        , 'Martin Roelsgaard Jakobsen'
        , 'Christian Kanstrup Holm'
        , 'Jens Christian Jensenius'
        , 'Steen Vang Petersen'
        , 'Steffen Thiel'
        , 'Trine Mogensen'
        , 'Søren Riis Paludan'
        , 'Tove Christensen'
        , 'Per Brøndsted Höllsberg'
        , 'Mart Ustav'
        , 'Andres Merits'
        , 'Malin Flodström-Tullberg'
        , 'Heikki Hyöty'
        , 'Veijo Hukkanen'
        , 'Tuomas Jartti'
        , 'Ilkka Julkunen'
        , 'Petri Susi'
        , 'Sisko Tauriainen'
        , 'Tytti Vuorinen'
        , 'Vincenzo Cerullo'
        , 'Andris Dišlers'
        , 'Kaspars Tārs'
        , 'Andris Zeltiņš'
        , 'Irina Sominska'
        , 'Louis Maes'
        , 'Koen Augustyns'
        , 'Kevin Ariën'
        , 'Karl Farrow'
        , 'Aya Takeoka'
        , 'Xavier Saelens'
        , 'Han Remaut'
        , 'Pascale QUATRESOOZ'
        , 'Danielle Perez Bercoff'
    ]

    for author in authors:
        # get the query for the API
        query = getQueryTerm(author)
        # call the search method and process the results
        search_results = search(query)
        # get the paper ids for the researcher and loop over it
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
        sleep(5)

    # Pretty print the first paper in full to observe its structure
    print(json.dumps(paper_summary, indent=2, separators=(',', ':')))
