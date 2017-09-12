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
**Tests:** Maybe it's enough, at least for now, to narrow down the search
field in the way that is already implemented and then apply the SQL cone search
function. I need to see how slow this is and then do more research.


# To Do List

1. Implement configuration file functionality
2. Implement cone search in SQL for extracting sources from tables within
circular region. Or implement alpha function in SQL for extraction from
square region on the sky. Since an SQL function needs to be defined either way,
the circular region is preferred.
3. Update database schema:
    - add Error and AssocSource tables
    - change Island table to rawIsland
    - change Source table to rawSource and remove cross-matching results
4. Implement source association
5. Implement image quality checks
6. Write testing scripts
7. Write documentation
