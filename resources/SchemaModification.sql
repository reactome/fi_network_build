-- Create a new class called "Interaction"
insert into `DataModel` values ("Interaction", "SchemaClass", "name", "Interaction", "STRING", 0);
insert into `DataModel` values ("Interaction", "SchemaClass", "super_classes", "Event", "SchemaClass", 0);
insert into `DataModel` values ("Interaction", "SchemaClass", "abstract", NULL, "SYMBOL", 0);

-- Create the Interaction table
DROP TABLE IF EXISTS `Interaction`;
CREATE TABLE `Interaction` (
	`DB_ID` int(10) unsigned NOT NULL,
	`interactionType` varchar(255),
	CONSTRAINT interaction_pk
		PRIMARY KEY (`DB_ID`)
) ENGINE = INNODB;

-- Create an attribute called "interactor"
insert into `DataModel` values ("interactor", "SchemaClassAttribute", "name", "interactor", "STRING", 0);
insert into `DataModel` values ("interactor", "SchemaClassAttribute", "min_cardinality", "2", "INTEGER", 0);
insert into `DataModel` values ("interactor", "SchemaClassAttribute", "allowed_classes", "PhysicalEntity", "SchemaClass", 0);
insert into `DataModel` values ("interactor", "SchemaClassAttribute", "max_carninality", NULL, "INTEGER", 0);
insert into `DataModel` values ("interactor", "SchemaClassAttribute", "multiple", "TRUE", "SYMBOL", 0);
insert into `DataModel` values ("interactor", "SchemaClassAttribute", "type", "db_instance_type", "STRING", 0);
insert into `DataModel` values ("interactor", "SchemaClassAttribute", "db_col_type", NULL, "STRING", 0);
insert into `DataModel` values ("Interaction:interactor", "SchemaClassAttribute", "name", "interactor", "STRING", 0);
insert into `DataModel` values ("Interaction:interactor", "SchemaClassAttribute", "allowed_classes", "PhysicalEntity", "SchemaClass", 0);
insert into `DataModel` values ("Interaction:interactor", "SchemaClassAttribute", "class", "Interaction", "SchemaClass", 0);
insert into `DataModel` values ("Interaction:interactor", "SchemaClassAttribute", "min_cardinality", "2", "INTEGER", 0);
insert into `DataModel` values ("Interaction:interactor", "SchemaClassAttribute", "max_carninality", NULL, "INTEGER", 0);
insert into `DataModel` values ("Interaction:interactor", "SchemaClassAttribute", "multiple", "TRUE", "SYMBOL", 0);
insert into `DataModel` values ("Interaction:interactor", "SchemaClassAttribute", "type", "db_instance_type", "STRING", 0);
insert into `DataModel` values ("Interaction:interactor", "SchemaClassAttribute", "db_col_type", NULL, "STRING", 0);
insert into `DataModel` values ("Interaction:interactor", "SchemaClassAttribute", "value_defines_instance", "all", "SYMBOL", 0);
insert into `DataModel` values ("Interaction:interactor", "SchemaClassAttribute", "category", "MANDATORY", "SYMBOL", 0);

-- Create an attribute called "interactionType"
insert into `DataModel` values ("interactionType", "SchemaClassAttribute", "name", "interactionType", "STRING", 0);
insert into `DataModel` values ("interactionType", "SchemaClassAttribute", "min_cardinality", "0", "INTEGER", 0);
insert into `DataModel` values ("interactionType", "SchemaClassAttribute", "max_carninality", "1", "INTEGER", 0);
insert into `DataModel` values ("interactionType", "SchemaClassAttribute", "multiple", "FALSE", "SYMBOL", 0);
insert into `DataModel` values ("interactionType", "SchemaClassAttribute", "type", "db_string_type", "STRING", 0);
insert into `DataModel` values ("interactionType", "SchemaClassAttribute", "db_col_type", NULL, "STRING", 0);
insert into `DataModel` values ("Interaction:interactionType", "SchemaClassAttribute", "name", "interactionType", "STRING", 0);
insert into `DataModel` values ("Interaction:interactionType", "SchemaClassAttribute", "class", "Interaction", "SchemaClass", 0);
insert into `DataModel` values ("Interaction:interactionType", "SchemaClassAttribute", "min_cardinality", "0", "INTEGER", 0);
insert into `DataModel` values ("Interaction:interactionType", "SchemaClassAttribute", "max_carninality", "1", "INTEGER", 0);
insert into `DataModel` values ("Interaction:interactionType", "SchemaClassAttribute", "multiple", "FALSE", "SYMBOL", 0);
insert into `DataModel` values ("Interaction:interactionType", "SchemaClassAttribute", "type", "db_string_type", "STRING", 0);
insert into `DataModel` values ("Interaction:interactionType", "SchemaClassAttribute", "db_col_type", NULL, "STRING", 0);
insert into `DataModel` values ("Interaction:interactionType", "SchemaClassAttribute", "category", "OPTIONAL", "SYMBOL", 0);

-- Create table for "interactor"
DROP TABLE IF EXISTS `Interaction_2_interactor`;
CREATE TABLE `Interaction_2_interactor` (
  `DB_ID` int(10) unsigned default NULL,
  `interactor_rank` int(10) unsigned default NULL,
  `interactor` int(10) unsigned default NULL,
  `interactor_class` varchar(64) default NULL,
  KEY `DB_ID` (`DB_ID`),
  KEY `interactor` (`interactor`)
) ENGINE=InnoDB;

-- Create a new attribute for DatabaseObject, dataSource. Its type should be ReferenceDatabase.
insert into `DataModel` values ("dataSource", "SchemaClassAttribute", "name", "dataSource", "STRING", 0);
insert into `DataModel` values ("dataSource", "SchemaClassAttribute", "min_cardinality", NULL, "INTEGER", 0);
insert into `DataModel` values ("dataSource", "SchemaClassAttribute", "max_cardinality", NULL, "INTEGER", 0);
insert into `DataModel` values ("dataSource", "SchemaClassAttribute", "multiple", "TRUE", "SYMBOL", 0);
insert into `DataModel` values ("dataSource", "SchemaClassAttribute", "type", "db_instance_type", "STRING", 0);
insert into `DataModel` values ("dataSource", "SchemaClassAttribute", "allowed_classes", "ReferenceDatabase", "SchemaClass", 0);
insert into `DataModel` values ("dataSource", "SchemaClassAttribute", "db_col_type", NULL, "STRING", 0);
insert into `DataModel` values ("DatabaseObject:dataSource", "SchemaClassAttribute", "name", "dataSource", "STRING", 0);
insert into `DataModel` values ("DatabaseObject:dataSource", "SchemaClassAttribute", "min_cardinality", NULL, "INTEGER", 0);
insert into `DataModel` values ("DatabaseObject:dataSource", "SchemaClassAttribute", "max_cardinality", NULL, "INTEGER", 0);
insert into `DataModel` values ("DatabaseObject:dataSource", "SchemaClassAttribute", "multiple", "TRUE", "SYMBOL", 0);
insert into `DataModel` values ("DatabaseObject:dataSource", "SchemaClassAttribute", "type", "db_instance_type", "STRING", 0);
insert into `DataModel` values ("DatabaseObject:dataSource", "SchemaClassAttribute", "allowed_classes", "ReferenceDatabase", "SchemaClass", 0);
insert into `DataModel` values ("DatabaseObject:dataSource", "SchemaClassAttribute", "db_col_type", NULL, "STRING", 0);
insert into `DataModel` values ("DatabaseObject:dataSource", "SchemaClassAttribute", "class", "DatabaseObject", "SchemaClass", 0);
insert into `DataModel` values ("DatabaseObject:dataSource", "SchemaClassAttribute", "category", "OPTIONAL", "SYMBOL", 0);

-- Need to create a new table for DatabaseObject:dataSource
DROP TABLE IF EXISTS `DatabaseObject_2_dataSource`;
CREATE TABLE `DatabaseObject_2_dataSource` (
  `DB_ID` int(10) unsigned default NULL,
  `dataSource_rank` int(10) unsigned default NULL,
  `dataSource` int(10) unsigned default NULL,
  `dataSource_class` varchar(64) default NULL,
  KEY `DB_ID` (`DB_ID`),
  KEY `dataSource` (`dataSource`)
) ENGINE=InnoDB;	

-- Increase the size of pages to 50. The original is too small.
ALTER TABLE LiteratureReference MODIFY pages varchar(50);

-- Increase the size of identifier to 50.
ALTER TABLE DatabaseIdentifier MODIFY identifier varchar(50);

-- create a new class called TargettedInteration for TF/Target or miRNA/Target interactions.
insert into `DataModel` values ("TargetedInteraction", "SchemaClass", "name", "TargetedInteraction", "STRING", 0);
insert into `DataModel` values ("TargetedInteraction", "SchemaClass", "super_classes", "Event", "SchemaClass", 0);
insert into `DataModel` values ("TargetedInteraction", "SchemaClass", "abstract", NULL, "SYMBOL", 0);

-- Create the TargetedInteraction table
DROP TABLE IF EXISTS `TargetedInteraction`;
CREATE TABLE `TargetedInteraction` (
	`DB_ID` int(10) unsigned NOT NULL,
	`factor` int(10) unsigned DEFAULT NULL,
	`factor_class` varchar(64) DEFAULT NULL,
	`target` int(10) unsigned DEFAULT NULL,
	`target_class` varchar(64) DEFAULT NULL,
	CONSTRAINT targettedInteraction_pk
		PRIMARY KEY (`DB_ID`)
) ENGINE = INNODB;

-- Create an attribute called "factor"
insert into `DataModel` values ("factor", "SchemaClassAttribute", "name", "factor", "STRING", 0);
insert into `DataModel` values ("factor", "SchemaClassAttribute", "min_cardinality", "0", "INTEGER", 0);
insert into `DataModel` values ("factor", "SchemaClassAttribute", "allowed_classes", "PhysicalEntity", "SchemaClass", 0);
insert into `DataModel` values ("factor", "SchemaClassAttribute", "max_carninality", "1", "INTEGER", 0);
insert into `DataModel` values ("factor", "SchemaClassAttribute", "multiple", "FALSE", "SYMBOL", 0);
insert into `DataModel` values ("factor", "SchemaClassAttribute", "type", "db_instance_type", "STRING", 0);
insert into `DataModel` values ("factor", "SchemaClassAttribute", "db_col_type", NULL, "STRING", 0);
insert into `DataModel` values ("TargetedInteraction:factor", "SchemaClassAttribute", "name", "factor", "STRING", 0);
insert into `DataModel` values ("TargetedInteraction:factor", "SchemaClassAttribute", "allowed_classes", "PhysicalEntity", "SchemaClass", 0);
insert into `DataModel` values ("TargetedInteraction:factor", "SchemaClassAttribute", "class", "TargetedInteraction", "SchemaClass", 0);
insert into `DataModel` values ("TargetedInteraction:factor", "SchemaClassAttribute", "min_cardinality", "0", "INTEGER", 0);
insert into `DataModel` values ("TargetedInteraction:factor", "SchemaClassAttribute", "max_carninality", "1", "INTEGER", 0);
insert into `DataModel` values ("TargetedInteraction:factor", "SchemaClassAttribute", "multiple", "FALSE", "SYMBOL", 0);
insert into `DataModel` values ("TargetedInteraction:factor", "SchemaClassAttribute", "type", "db_instance_type", "STRING", 0);
insert into `DataModel` values ("TargetedInteraction:factor", "SchemaClassAttribute", "db_col_type", NULL, "STRING", 0);
insert into `DataModel` values ("TargetedInteraction:factor", "SchemaClassAttribute", "value_defines_instance", "all", "SYMBOL", 0);
insert into `DataModel` values ("TargetedInteraction:factor", "SchemaClassAttribute", "category", "MANDATORY", "SYMBOL", 0);

-- Create an attribute called "target"
insert into `DataModel` values ("target", "SchemaClassAttribute", "name", "target", "STRING", 0);
insert into `DataModel` values ("target", "SchemaClassAttribute", "min_cardinality", "0", "INTEGER", 0);
insert into `DataModel` values ("target", "SchemaClassAttribute", "allowed_classes", "PhysicalEntity", "SchemaClass", 0);
insert into `DataModel` values ("target", "SchemaClassAttribute", "max_carninality", "1", "INTEGER", 0);
insert into `DataModel` values ("target", "SchemaClassAttribute", "multiple", "FALSE", "SYMBOL", 0);
insert into `DataModel` values ("target", "SchemaClassAttribute", "type", "db_instance_type", "STRING", 0);
insert into `DataModel` values ("target", "SchemaClassAttribute", "db_col_type", NULL, "STRING", 0);
insert into `DataModel` values ("TargetedInteraction:target", "SchemaClassAttribute", "name", "target", "STRING", 0);
insert into `DataModel` values ("TargetedInteraction:target", "SchemaClassAttribute", "allowed_classes", "PhysicalEntity", "SchemaClass", 0);
insert into `DataModel` values ("TargetedInteraction:target", "SchemaClassAttribute", "class", "TargetedInteraction", "SchemaClass", 0);
insert into `DataModel` values ("TargetedInteraction:target", "SchemaClassAttribute", "min_cardinality", "0", "INTEGER", 0);
insert into `DataModel` values ("TargetedInteraction:target", "SchemaClassAttribute", "max_carninality", "1", "INTEGER", 0);
insert into `DataModel` values ("TargetedInteraction:target", "SchemaClassAttribute", "multiple", "FALSE", "SYMBOL", 0);
insert into `DataModel` values ("TargetedInteraction:target", "SchemaClassAttribute", "type", "db_instance_type", "STRING", 0);
insert into `DataModel` values ("TargetedInteraction:target", "SchemaClassAttribute", "db_col_type", NULL, "STRING", 0);
insert into `DataModel` values ("TargetedInteraction:target", "SchemaClassAttribute", "value_defines_instance", "all", "SYMBOL", 0);
insert into `DataModel` values ("TargetedInteraction:target", "SchemaClassAttribute", "category", "MANDATORY", "SYMBOL", 0);
