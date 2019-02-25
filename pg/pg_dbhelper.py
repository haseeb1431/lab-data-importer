import psycopg2
from pg.config import config
import requests
import json


def insertData(sqlStatement):
    if sqlStatement is None:
        return False

    conn = None
    creds = config()

    try:
        # establish the connection with the server
        conn = psycopg2.connect(**creds)

        # create a cursor
        cur = conn.cursor()
        # execute the sql statement on the server
        cur.execute(sqlStatement)
        # cur.executeMany(sqlStatement,data dictionary)
        conn.commit()

        return True

    except (Exception, psycopg2.DatabaseError) as error:
        print(error)

    finally:
        if conn is not None:
            # print("Closing connection")
            conn.close()

    return False


def getRowSql(row):
    sql = ""
    if row is not None:
        sql = " ('{0}','{1}','{2}','{3}','{4}','{5}','{6}','{7}')".format(row["FullJournalName"],
                                                                          row["PubTypeList"], \
                                                                          row["LastAuthor"],
                                                                          row["AuthorList"],
                                                                          row["PubDate"], \
                                                                          row["Title"],
                                                                          row["Id"],
                                                                          row["PI"])

    return sql


def getBooks():
    # hardcoding the url instead of config driven
    url = "https://api.nytimes.com/svc/books/v3/lists/names.json?api-key=9324744857a342b9a546a055cbd8a578"

    # Send the get request
    myResponse = requests.get(url)

    # For successful API call, response code will be 200 (OK)
    if (myResponse.ok):

        # parse the json content
        jData = json.loads(myResponse.content)

        print("We found {0} rows \n".format(len(jData["results"])))

        return jData["results"]

    #        insertSql= 'INSERT INTO public.nybooks( display_name, list_name, list_name_encoded, newest_published_date, oldest_published_date, updated)	VALUES ';
    #        valuesSQL = ""
    #        #valuesSQL = getBookSql(jData["results"][0])
    #        #best is to use bulk write
    #
    #        for row in jData["results"]:
    #            valuesSQL+=  ",{0}".format(getBookSql(row))
    #
    #        insertSql += valuesSQL[1:len(valuesSQL)] # remove the first comma
    #        #print(insertSql)
    #
    #        InsertData(insertSql)
    # print("Data is uploaded successfully")
    else:
        # If response code is not ok (200), print the resulting http error code with description
        myResponse.raise_for_status()


def generateSQL(paper_summary):
    if paper_summary is not None:

        insertSql = 'INSERT INTO public."Jsummary"( "FullJournalName", "PubTypeList", "LastAuthor", "AuthorList", "PubDate", "Title", "PMID","PIAuthor")	VALUES ';
        valuesSQL = ""

        for row in paper_summary:
            valuesSQL += ",{0}".format(getRowSql(row))

        insertSql += valuesSQL[1:len(valuesSQL)]  # remove the first comma
        # print(insertSql)

        return insertSql

    return None
