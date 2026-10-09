# search/v422: frozen backup of engine v4.22

A copy of `search/` at commit ad32f75 (engine v4.22), kept as a backup when engine v9 replaced it on 2026-10-08.
Not selectable by main. Imports point at this package (`search.v422.*`) so it stays importable for reference.
`test_engine.py` is its test file (not collected by `python -m unittest discover -s tests`).
History: `DEVELOPMENT.md` and `README.md` here.
Its settings (opening, domains, schedule, merges, sharing, scout, riders, climb, pick) were removed from
`common/tuning.py` for v9, so it no longer runs as is: check out ad32f75 to run it.
It reads the portfolio files from search/portfolio/ (the same ones v9 reads), not a copy.
