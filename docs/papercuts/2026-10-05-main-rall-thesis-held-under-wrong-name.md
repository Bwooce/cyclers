# Rall's Sc.D. thesis was held for months as "hollister-rall-1970-periodic-orbits-NASA-CR.pdf"

- Date: 2026-10-05
- Agent: main
- Seen before: yes (processed/2026-10-05-main-corpus-index-abbreviated-filenames.md; the Bruno/Brjuno miss in the wanted-list review)

What happened: the wanted list asked the owner for Rall 1969 (PhD, MIT); the lead fetched a 219-page image-only scan from DSpace and queued it for OCR. The owner's NTRS 19700017824 is byte-identical to the held "hollister-rall-1970 NASA CR", whose title page is Rall's thesis (MIT MSL report TE-34) with a full text layer. check_wanted_vs_corpus.py missed it: the filename's first author is the supervisor.
Workaround: md5 match on the owner's upload; corpus-file-opus told to correct index, notes and wanted list.
Suggested fix: index rows carry title-page author and title (not the filer's guess) and report/thesis numbers; the checker matches on those fields too.
