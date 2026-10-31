import sys, fitz
doc = fitz.open(sys.argv[1])
open(sys.argv[2], "w", encoding="utf-8").write("\n".join(f"=== page {i+1}\n" + p.get_text() for i, p in enumerate(doc)))
print(len(doc))
