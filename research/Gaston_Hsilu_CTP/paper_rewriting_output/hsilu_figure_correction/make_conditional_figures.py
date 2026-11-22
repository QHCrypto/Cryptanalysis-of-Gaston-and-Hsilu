#!/usr/bin/env python3
"""Transcribe the supplied SVG paths, validate degree propagation, and draw TikZ."""
from pathlib import Path
import collections
import hashlib
import json
import re
import sys
import xml.etree.ElementTree as ET

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[1]
ROOT = PAPER.parent
sys.path.insert(0, str(ROOT/'代数攻击/hsilu_2r'))
from verify_path_algebra import PARAM

PALETTE = {'#ffd84d': 'aff', '#b8b8b8': 'fixed', '#d73a49': 'zero',
           '#2da44e': 'one', '#8ec5ff': 'selected', '#ffffff': 'other',
           'white': 'other'}
COLORS = {'aff':'F4D77B', 'fixed':'D7DBDF', 'zero':'E89494',
          'one':'8CC9A2', 'selected':'8CB9E3', 'other':'FFFFFF'}

def read_svg(rounds):
    src = ROOT/('代数攻击/hsilu_2r/hsilu_no_start.svg' if rounds == 2
                else '代数攻击/Hsilu3r_new/hsilu_3round_best_attack_1_5.1.svg')
    names = sum(([f'X{r}', f'Y{r}'] for r in range(1, rounds+1)), []) + [f'X{rounds+1}']
    states = {name:[[None]*64 for _ in range(5)] for name in names}
    for e in ET.parse(src).getroot().findall('{http://www.w3.org/2000/svg}rect'):
        title = e.find('{http://www.w3.org/2000/svg}title')
        match = re.fullmatch(r'([XY]\d+)\[(\d+)\]\[(\d+)\]', title.text or '') if title is not None else None
        if match:
            name, row, col = match.group(1), int(match.group(2)), int(match.group(3))
        elif e.get('width') == '7' and e.get('height') == '7':
            x,y = int(e.get('x')), int(e.get('y'))
            stage, remainder = divmod(y-78, 78)
            name, row, col = names[stage], remainder//8, (x-86)//8
        else:
            continue
        assert states[name][row][col] is None
        states[name][row][col] = PALETTE[e.get('fill')]
    assert all(c is not None for a in states.values() for row in a for c in row)
    return src, states

def add(x,y):
    if x[1] is not None and y[1] is not None:
        return (0, x[1]^y[1])
    return (max(x[0],y[0]), None)

def mul(x,y):
    if x[1] == 0 or y[1] == 0: return (0,0)
    if x[1] == 1: return y
    if y[1] == 1: return x
    return (min(2,x[0]+y[0]), None)

def linear(a):
    x=list(a)
    for k in range(5,18):
        x.append([add(x[PARAM['n1'][k]][j],x[PARAM['n2'][k]][(j+PARAM['m'][k])%64]) for j in range(64)])
    return [[x[13+i][(j+PARAM['s'][i])%64] for j in range(64)] for i in range(5)]

def check(rounds, states):
    active=[j for j,c in enumerate(states['X1'][0]) if c=='aff']
    x=[[(int(i==0 and j in active),None) for j in range(64)] for i in range(5)]
    verified=0
    conditions={}
    def compare(name, a):
        nonlocal verified
        for i,row in enumerate(a):
            for j,(degree,known) in enumerate(row):
                code=states[name][i][j]
                expect=0 if code in ('fixed','zero','one') else 1 if code in ('aff','selected') else 2
                assert degree == expect, (name,i,j,code,degree)
                verified+=1
    for r,rc in enumerate([0x69,0x5a,0x4b][-rounds:],1):
        conditions[f'X{r}']=[]
        for i in range(5):
            for j in range(64):
                c=states[f'X{r}'][i][j]
                if c in ('zero','one'):
                    assert x[i][j][0] == 0
                    value=int(c=='one')
                    x[i][j]=(0,value)
                    conditions[f'X{r}'].append([i,63-j,value])
        compare(f'X{r}',x)
        for j in range(64):
            if rc>>(63-j)&1: x[0][j]=add(x[0][j],(0,1))
        y=[[add(x[i][j],mul(add(x[(i+1)%5][j],(0,1)),x[(i+2)%5][j])) for j in range(64)] for i in range(5)]
        compare(f'Y{r}',y)
        x=linear(y)
    compare(f'X{rounds+1}',x)
    return {'verified_cells':verified,'conditions_native':conditions,
            'active_native':sorted(63-j for j in active),
            'selected_native':sorted(63-j for j,c in enumerate(states[f'X{rounds+1}'][0]) if c=='selected')}

def draw(rounds, states, report):
    lines=[r'\begingroup']
    lines += [rf'\definecolor{{hcf{k}}}{{HTML}}{{{v}}}' for k,v in COLORS.items()]
    lines += [r'\resizebox{\textwidth}{!}{%',
              r'\begin{tikzpicture}[x=1.9mm,y=-1.9mm,>=Stealth,line width=.25pt,every node/.style={font=\scriptsize}]']
    def text(x,y,s,anchor='west',font=r'\scriptsize'):
        lines.append(rf'\node[anchor={anchor},font={font}] at ({x},{y}) {{{s}}};')
    lines.append(r'\path (-5,-1) -- (85,-1);')
    legend=[(0,0,'aff','affine'),(32,0,'fixed','fixed in the structure'),
            (0,3,'selected','selected affine output'),(32,3,'other','other coordinates'),
            (0,6,'zero','condition $=0$'),(32,6,'one','condition $=1$')]
    for x,y,c,label in legend:
        lines.append(rf'\filldraw[fill=hcf{c},draw=black!55] ({x},{y}) rectangle ++(1.3,1.3);')
        text(x+2,y+.65,label)
    names=list(states)
    for stage,name in enumerate(names):
        y=12+11*stage
        label=rf'$X^{{({stage//2})}}$' if stage%2==0 else rf'$U^{{({stage//2+1})}}$'
        text(-1,y+2.5,label,'east')
        for j in [0,8,16,24,32,40,48,56,63]: text(10+j+.5,y-.6,str(j),'south',r'\tiny')
        for i in range(5):
            text(9.3,y+i+.5,str(i),'east',r'\tiny')
            for j in range(64):
                c=states[name][i][63-j]
                lines.append(rf'\fill[hcf{c}] ({10+j},{y+i}) rectangle ++(1,1);')
                if c in ('zero','one'):
                    text(10+j+.5,y+i+.51,str(int(c=='one')),'center',r'\fontsize{4}{4}\selectfont')
        for j in range(65): lines.append(rf'\draw[black!45] ({10+j},{y}) -- ({10+j},{y+5});')
        for i in range(6): lines.append(rf'\draw[black!45] (10,{y+i}) -- (74,{y+i});')
        if stage < len(names)-1:
            lines.append(rf'\draw[->,line width=.6pt] (40,{y+5.5}) -- (40,{y+8.5});')
            lab=rf'$C_{{{12-rounds+stage//2}}},\ \chi$' if stage%2==0 else r'$L_\Hh$'
            text(42,y+7,lab)
    footer=rf'{len(report["active_native"])} active variables; '
    footer+='; '.join(f'{len(v)} condition{"s" if len(v)!=1 else ""} on $X^{{({int(k[1:])-1})}}$' for k,v in report['conditions_native'].items())
    text(37,12+11*(len(names)-1)+9,footer,'center')
    lines += [r'\end{tikzpicture}}',r'\endgroup']
    code='\n'.join(lines)+'\n'
    name=f'hsilu{rounds}_preimage'
    for folder in ['manuscript','translation_comparison']:
        out=PAPER/folder/'figures';out.mkdir(exist_ok=True)
        (out/f'{name}.tex').write_text(code)
    wrapper=r'''\documentclass[border=4pt]{standalone}
\usepackage{amsmath,amssymb,graphicx,tikz}
\usetikzlibrary{arrows.meta}
\newcommand{\Hh}{\mathsf H}
\begin{document}\setlength{\textwidth}{170mm}
\input{'''+name+r'''.tex}
\end{document}
'''
    (PAPER/'manuscript/figures'/f'{name}_standalone.tex').write_text(wrapper)

def main():
    result={}
    for rounds in (2,3):
        src,states=read_svg(rounds)
        report=check(rounds,states)
        report.update(source=str(src), source_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),
                      convention='source MSB=0; manuscript native LSB=0; native=63-source', status='PASS')
        draw(rounds,states,report)
        result[str(rounds)]=report
    (HERE/'figure_validation.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:{x:v[x] for x in ['verified_cells','active_native','selected_native','status']} for k,v in result.items()},indent=2))

if __name__=='__main__': main()
