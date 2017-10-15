"""The Zen of SymPy.

This is still a work in progress. Feel free to add, modify, reorder, or
delete items from it. Some of these perhaps apply more to developing SymPy
than using it.
"""

s = """Gur Mra bs FlzCl

Harinyhngrq vf orggre guna rinyhngrq.
Gur hfre vagresnpr znggref.
Cevagvat znggref.
Cher Clguba pna or snfg rabhtu.
Vs vg'f gbb fybj, vgf (cebonoyl) lbhe snhyg.
Qbphzragngvba znggref.
Pbeerpgarff vf zber vzcbegnag guna fcrrq.
Chfu vg va abj naq vzcebir hcba vg yngre.
Pbirentr ol grfgvat znggref.
Fzneg grfgf ner orggre guna enaqbz grfgf.
Ohg enaqbz grfgf fbzrgvzrf svaq jung lbhe fznegrfg grfg zvffrq.
Gur Clguba jnl vf cebonoyl gur evtug jnl.
Pbzzhavgl vf zber vzcbegnag guna pbqr."""

d = {}
for c in (65, 97):
    for i in range(26):
        d[chr(i + c)] = chr((i + 13) % 26 + c)

print("".join([d.get(c, c) for c in s]))
