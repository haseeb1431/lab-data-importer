--select * from public."Jsummary"


select "PIAuthoer, count(PIAuthor)
From public."Jsummary"
Group by PIAuthoer



Select 
	js."FullJournalName",
	jif."JournalImpact",
	js."Title",
	js."LastAuthor"
From
public."JIF2018" jif
inner join public."Jsummary" js on jif."FullJournalTitle"=js."FullJournalName"