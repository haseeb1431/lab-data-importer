CREATE TABLE public."JIF2018"
(
    "Rank" serial NOT NULL,
    "FullJournalTitle" character varying(500) NOT NULL,
    "TotalCites" bigint,
    "JournalImpact" double precision NOT NULL,
    "Eigenfactor score" double precision,
    PRIMARY KEY ("Rank")
)
WITH (
    OIDS = FALSE
);

ALTER TABLE public."JIF2018"
    OWNER to postgres;