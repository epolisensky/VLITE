# Ongoing Considerations

#### SQLite or Client/Server RDBMS
**Currently using:** SQLite

**SQLite benefits:**
- Comes with Python
- Single file on the system; don't have to setup database server
- Fast

**Client/server RDBMS benefits:**
- Allows concurrent writes -- not an issue for us
- Better for really large databases -- might become an issue for us
- Allows our database to be accessed by people over a network (i.e. like SDSS)

**Tests:** None for now

#### Raw SQL vs. ORM (SQLAlchemy)
**Currently using:** raw SQL

**Raw SQL benefits:**
- Possibly faster performance
- Already coded

**SQLAlchemy benefits:**
- Makes our code portable by avoiding having to re-write SQL statements if we do
decide to switch from SQLite to a client/server RDBMS
- Reduces complexity?
- Can still use raw SQL statements if I can't figure out how to do it
better/faster in SQLAlchemy

**Tests:** Just going to implement both ways and test which is faster and more
convenient.

#### Spatial indexing for faster database searches
**Currently using:** nothing

**Benefits:** Although not yet implemented (see to do list item #2.), the SQL
queries for extracting sources within a region on the sky can be pretty costly
(at least based on what I've read). To avoid having to scan large tables, many
astronomical databases employ some form of spatial indexing to divide sources
into sky regions to narrow down the search range.

**Disadvantages:** Time it will take me to fully understand this problem and
figure out how to implement it.

**Tests:** Not sure.


# To Do List

- [x] Implement configuration file functionality
- [x] Implement cone search in SQL for extracting sources from tables within
circular region. 
- [x] Update database schema
- [ ] Implement source association
- [ ] Update catalog cross-matching
- [ ] Implement image quality checks
- [ ] Write testing scripts
- [ ] Write documentation