create table categories (
  id integer not null,
  catname varchar(40) not null,
  primary key(id)
);

create table authors (
  id integer not null,
  autname varchar(100) not null,
  primary key(id)
);

create table quotes (
  id      integer not null,
  cid     integer not null,
  aid     integer,
  quotext varchar(250) not null,
  primary key(id)
);


