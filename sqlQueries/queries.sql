--select * from public."Jsummary"


Select 
	"PIAuthor", count("PIAuthor") as cnt
FROM public."Jsummary"
Group by "PIAuthor"
Order by 2 desc


-- Merge
Select 
	js."FullJournalName",
	jif."JournalImpact",
	js."Title",
	js."LastAuthor"
From
public."JIF2018" jif
inner join public."Jsummary" js on jif."FullJournalTitle"=js."FullJournalName"

-- sum 
Select 
	
	sum(jif."JournalImpact") as totalIF,
	Count(jif."JournalImpact") as totalPUblications,
	js."LastAuthor"
From
public."JIF2018" jif
inner join public."Jsummary" js on jif."FullJournalTitle"=js."FullJournalName"
Group by js."LastAuthor"