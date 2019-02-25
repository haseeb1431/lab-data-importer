import Bio

def parse_Journal(journal):
    """            <Journal>
                <ISSN IssnType="Electronic">1573-7276</ISSN>
                <JournalIssue CitedMedium="Internet">
                    <PubDate>
                        <Year>2019</Year>
                        <Month>Feb</Month>
                        <Day>04</Day>
                    </PubDate>
                </JournalIssue>
                <Title>Clinical &amp; experimental metastasis</Title>
                <ISOAbbreviation>Clin. Exp. Metastasis</ISOAbbreviation>
            </Journal>"""



def parse_keywords(keyword_list):
    """Parse keywords from article """
    keywords = list()
    if keyword_list is not None:
        if type(keyword_list) is list:
            keyword_list = keyword_list

        for k in keyword_list:
            if str(k) is not None:
                keywords.append(str(k))
        keywords = '; '.join(keywords)
    else:
        keywords = ''
    return keywords



def parse_publication_type(PublicationTypeList):

    publication_type = list()
    if PublicationTypeList is not None:
        # if type(PublicationTypeList) is Bio.Entrez.Parser.ListElement:
        #     PublicationTypeList = PublicationTypeList[0]

        for k in PublicationTypeList:
            if str(k) is not None:
                publication_type.append(str(k))
        publication_type = '; '.join(publication_type)
    else:
        publication_type = ''

    return publication_type