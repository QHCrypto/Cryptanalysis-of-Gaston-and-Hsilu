#!/usr/bin/env python3
"""Derive restricted-kernel matrices from the official round's linear map."""
import ctypes
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile
from audit_kernel_structure import maps, rol, EAST, WEST
from verify_lossy_model import chi, column

ROOT=Path(__file__).resolve().parents[1]
REFERENCE=ROOT.parent/'gaston_ref/gaston/implementation/reference'
class State(ctypes.Structure):
    _fields_=[('x',ctypes.c_uint64*5)]

inverse={chi(x):x for x in range(32)}
tables=[]
with tempfile.TemporaryDirectory(prefix='gaston-support-reference-') as temp:
    path=Path(temp)/'reference.so'
    subprocess.run(['cc','-O2','-shared','-fPIC',str(REFERENCE/'gaston.c'),'-o',str(path)],check=True)
    lib=ctypes.CDLL(str(path));lib.gaston_round.argtypes=[ctypes.POINTER(State),ctypes.c_uint64]
    lib.gaston_round.restype=None
    def linear(words):
        a=State((ctypes.c_uint64*5)(*words));lib.gaston_round(ctypes.byref(a),0)
        columns=[inverse[column(a.x,z)] for z in range(64)]
        return [sum(((v>>i)&1)<<z for z,v in enumerate(columns)) for i in range(5)]
    for side in ('left','right'):
        table=[]
        for z in range(64):
            row=[]
            for i in range(5):
                unit=[0]*5;unit[i]=1<<z
                b=[rol(w,-EAST[j]-WEST[j]) for j,w in enumerate(unit)] if side=='left' else unit
                expected=[rol(w,EAST[j]+WEST[j]) for j,w in enumerate(b)]
                actual=linear(b)
                effect=rol(actual[0]^expected[0],-23)
                assert all(actual[j]^expected[j]==rol(effect,23+WEST[j]) for j in range(5))
                theta=[rol(w,EAST[j]) for j,w in enumerate(b)]
                assert effect==maps(theta)[3]
                row.append(effect)
            table.append(row)
        tables.append(table)
output=ROOT/'results/support_reference_matrices_20261024.txt'
output.write_text('\n'.join(' '.join(f'{v:016x}' for v in row) for table in tables for row in table)+'\n')
report={'status':'PASS','official_reference_basis_cases':640,'matrix_file':str(output),
        'matrix_sha256':hashlib.sha256(output.read_bytes()).hexdigest(),
        'reference_sha256':hashlib.sha256((REFERENCE/'gaston.c').read_bytes()).hexdigest(),
        'scope':'Left tests E(rho_west^-1(a1)); right tests E(rho_east(b1)); derived via the official full linear layer'}
(ROOT/'results/support_reference_verification_20261024.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
