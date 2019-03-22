DROP TABLE IF EXISTS Tests;
DROP TABLE IF EXISTS Isolates;
DROP TABLE IF EXISTS Microbes;
DROP TABLE IF EXISTS Samples;

CREATE TABLE Samples
(
	ID VARCHAR(15) NOT NULL PRIMARY KEY,
	SellByDate DATE,
	Month INT,
	Year INT,
	State VARCHAR(15),
	MeatType VARCHAR(15),
	Organic BOOLEAN,
	Agency VARCHAR(15),
	Source VARCHAR(15),
	HostSpecies VARCHAR(15)
);

CREATE TABLE Microbes
(
	ID INT NOT NULL PRIMARY KEY,
	Genus VARCHAR(5),
	GenusName VARCHAR(25),
	Species VARCHAR(25),
	Serotype VARCHAR(25),
	AntigenicFormula VARCHAR(25)
);

CREATE TABLE Isolates
(
	ID VARCHAR(15) NOT NULL PRIMARY KEY,
	Plate VARCHAR(15),
	SampleID VARCHAR(15),
	MicrobeID INT,
	CONSTRAINT Isolates_Samples_ID_fk FOREIGN KEY(SampleID) REFERENCES Samples (ID),
	CONSTRAINT Isolates_Microbes_ID_fk FOREIGN KEY(MicrobeID) REFERENCES Microbes (ID)
);

CREATE TABLE Tests
(
	IsolateID VARCHAR(15),
	Drug VARCHER(10),
	MIC DEC(8,4),
	SIR CHAR(1),
	Sign CHAR(2),
	PRIMARY KEY (IsolateID, Drug),
	CONSTRAINT Tests_Isolates_ID_fk FOREIGN KEY(IsolateID) REFERENCES Isolates (ID)
);


























