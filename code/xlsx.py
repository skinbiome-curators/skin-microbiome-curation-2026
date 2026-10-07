import zipfile, re
from xml.etree.ElementTree import iterparse

NS = '{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'

def shared_strings(z):
    out = []
    if 'xl/sharedStrings.xml' not in z.namelist():
        return out
    cur = []
    for ev, el in iterparse(z.open('xl/sharedStrings.xml'), events=('end',)):
        if el.tag == NS + 't':
            cur.append(el.text or '')
        elif el.tag == NS + 'si':
            out.append(''.join(cur)); cur = []; el.clear()
    return out

def col_idx(ref):
    m = re.match(r'([A-Z]+)', ref)
    n = 0
    for ch in m.group(1):
        n = n * 26 + (ord(ch) - 64)
    return n - 1

def rows(path, sheet='xl/worksheets/sheet1.xml', limit=None):
    z = zipfile.ZipFile(path)
    ss = shared_strings(z)
    cnt = 0
    for ev, el in iterparse(z.open(sheet), events=('end',)):
        if el.tag != NS + 'row':
            continue
        vals = {}
        for c in el.findall(NS + 'c'):
            t = c.get('t'); v = c.find(NS + 'v')
            if t == 'inlineStr':
                isel = c.find(NS + 'is')
                txt = ''.join(x.text or '' for x in isel.iter(NS + 't')) if isel is not None else ''
            elif v is None:
                txt = ''
            elif t == 's':
                txt = ss[int(v.text)]
            else:
                txt = v.text or ''
            vals[col_idx(c.get('r'))] = txt
        n = (max(vals) + 1) if vals else 0
        yield [vals.get(i, '') for i in range(n)]
        el.clear()
        cnt += 1
        if limit and cnt >= limit:
            return
