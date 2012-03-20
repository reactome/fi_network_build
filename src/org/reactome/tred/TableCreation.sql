create table `Source` ( 
`source_id` int(10) not null auto_increment,
`source_name` varchar(50) not null,
`source_version` varchar(10) default null,
`source_date` date default null,
`description` varchar(255) default null,
`base_uri` varchar(255) default null,
primary key(`source_id`)
)engine=MyISAM;

