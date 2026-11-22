#!/usr/bin/env python3
"""HNW26 full/zero local-column estimates for a saved best full-DATF endpoint."""
import argparse
import concurrent.futures
import csv
import hashlib
import html
import json
import subprocess
import time
from fractions import Fraction
from pathlib import Path

from fixed_tail_search import S,column
from hnw26_compare import ROOT,run_job,save,fraction
from hnw26_audit import reference_zero
from datf_linear_two_report import number,bits,states


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,required=True,help='Matching result directory containing endpoint_check/results.json and sampling.csv')
    args=p.parse_args();source=args.source.resolve();out=source/'hnw26';out.mkdir(exist_ok=True)
    path=source/'endpoint_check/results.json'
    selected=json.loads(path.read_text())['best_full_DATF_within_selection']
    b=selected['candidate'];eid=selected['endpoint_id'];datf=selected['result']['correlation']
    samples=list(csv.DictReader((source/'sampling.csv').open()))
    sample=next(s for s in samples if s['input_words'].split(';')==b['input_difference']
                and s['mask_words'].split(';')==b['output_mask']
                and int(s['rounds'])==b['R'] and int(s['start'])==b['start_round_index'])
    job={'id':eid,'candidate':b}
    sources=[path,source/'sampling.csv',ROOT/'src/hnw_endpoint.cpp',ROOT/'build/hnw_endpoint',
             ROOT/'include/gaston.hpp',ROOT/'include/estimator.hpp',Path(__file__)]
    config={'endpoint_id':eid,'candidate':b,'modes':['full','zero'],'rd':1,'rl':0,
            'first_L_omitted':True,'reference_constants_retained':True,
            'method':'HNW26 local-column approximation with exact first-round DDT enumeration',
            'full':'1024 joint state/difference Fourier coordinates per column',
            'zero':'32 zero-state-mask difference coordinates per column',
            'output_readout':'Product over active output columns within each first-DDT branch, then signed probability-weighted sum.',
            'scope':'One fixed complete endpoint; original 1+1+2 path weights are not multiplied into HNW again.',
            'sources':[{'path':str(f),'sha256':hashlib.sha256(f.read_bytes()).hexdigest()} for f in sources]}
    cfg=out/'configuration.json'
    if cfg.exists():assert json.loads(cfg.read_text())==config,'Configuration changed; choose a new source directory.'
    else:save(cfg,config)
    began=time.perf_counter()
    backend=json.loads(subprocess.check_output([str(ROOT/'build/hnw_endpoint'),'--selftest'],text=True))
    assert backend['status']=='PASS';save(out/'backend_validation.json',backend)
    estimates={}
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:
        futures={ex.submit(run_job,job,m,out):m for m in ('full','zero')}
        for f in concurrent.futures.as_completed(futures):
            mode=futures[f];estimates[mode]=f.result()
            print(mode,estimates[mode]['correlation'],estimates[mode]['absorbed_histogram'],flush=True)
    delta=tuple(int(w,16) for w in b['input_difference'])
    ddt=[[sum(S[x]^S[x^a]==v for x in range(32)) for v in range(32)] for a in range(32)]
    branch_checks=0
    artifacts={}
    for mode,r in estimates.items():
        raw=json.loads(Path(r['artifact']).read_text());artifacts[mode]=raw;seen=set()
        for br in raw['branches']:
            d=tuple(int(w,16) for w in br['first_output_difference'].split(';'))
            assert d not in seen;seen.add(d);prob=Fraction(1)
            for j in range(64):prob*=Fraction(ddt[column(delta,j)][column(d,j)],32)
            assert prob>0 and prob==fraction(br['probability']);branch_checks+=1
        assert fraction(r['bias'])*2==fraction(r['correlation'])
    oracle=reference_zero(b)
    assert oracle==fraction(estimates['zero']['correlation'])
    validation={'status':'PASS','DDT_branch_probabilities_independently_checked':branch_checks,
                'probability_mass_and_signed_contribution_sums':'PASS',
                'independent_Python_Fraction_zero_model':str(oracle),
                'full_and_zero_agree':fraction(estimates['full']['correlation'])==fraction(estimates['zero']['correlation']),
                'backend_selftest':backend,'unknown_or_truncated_results':False}
    save(out/'validation.json',validation)
    result={'status':'completed','endpoint_id':eid,'candidate':b,'DATF':datf,'HNW':estimates,
            'sampling':sample,'first_L_omitted':True,'rd':1,'rl':0,'elapsed_seconds':time.perf_counter()-began,
            'interpretation':'HNW model estimates, not exact full permutation correlations.'}
    save(out/'results.json',result)
    branch_rows=[]
    zero_lookup={x['first_output_difference']:x for x in artifacts['zero']['branches']}
    for i,x in enumerate(artifacts['full']['branches'],1):
        z=zero_lookup[x['first_output_difference']]
        branch_rows.append(f'<tr><td>{i}</td><td>{html.escape(bits(x["first_output_difference"].split(";")))}</td>'
                          f'<td>{number(x["probability"])}</td><td>{number(x["conditional_model_correlation"])}</td>'
                          f'<td>{number(z["conditional_model_correlation"])}</td><td>{x["absorbed_by_round"] or "未完全退化"}</td>'
                          f'<td>{number(x["contribution"])}</td></tr>')
    full=estimates['full'];zero=estimates['zero']
    all_branch_zero=full['nonzero_branch_correlations']==zero['nonzero_branch_correlations']==0
    cause=('两版的每个首轮分支最终均为模型零；不是正负非零分支求和抵消。'
           if all_branch_zero else '最终值为各首轮分支的条件模型相关度按精确 DDT 概率带符号求和。')
    model=float(fraction(full['correlation']));lo=float(sample['ci95_low']);hi=float(sample['ci95_high'])
    conclusion=('HNW 结果落在已有采样的 95% 区间内。' if lo<=model<=hi else
                'HNW 结果落在已有采样的 95% 区间之外；本端点上逐列近似没有准确反映采样观测到的完整相关度。')
    if all_branch_zero and lo>0:
        conclusion+=' 这里的 0 是近似模型内的零，不能据此否定该端点的非零现象。'
    page=f'''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Gaston 四轮新端点 · HNW26 复算</title><style>
body{{margin:0;background:#f3f6f4;color:#19392f;font:16px/1.8 system-ui,sans-serif}}main{{max-width:1150px;margin:auto;padding:32px 24px 72px}}h1{{font-size:32px;line-height:1.4}}h2{{font-size:23px}}section{{background:white;border:1px solid #d5e3dd;border-radius:10px;padding:24px;margin:24px 0}}a{{color:#087868}}table{{width:100%;border-collapse:collapse;font-size:14px}}th,td{{text-align:left;padding:12px 9px;border-bottom:1px solid #d5e3dd}}th{{background:#ecf3ef}}.scroll{{overflow:auto}}.callout{{background:#fff3df;border-left:4px solid #b67a1e;padding:16px 20px}}.states{{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:16px}}pre{{font:13px/1.8 monospace;overflow:auto;background:#eef4f0;padding:14px}}.muted{{color:#5b7267}}@media(max-width:750px){{main{{padding:20px 14px}}.states{{grid-template-columns:1fr}}section{{padding:18px}}}}</style></head><body><main>
<h1>Gaston 四轮新端点：HNW26 相关度复算</h1><p>固定刚才找到的输入差分与输出掩码；仍然只省略完整区间首个线性层。原路径来自 1+1+2，HNW 按此前的 r<sub>d</sub>=1、r<sub>l</sub>=0 设置重评整个四轮区间，即 1+3+0。</p>
<section><h2>固定端点与模型对照</h2><p>Δ = {{{html.escape(bits(b['input_difference']))}}}；Λ = {{{html.escape(bits(b['output_mask']))}}}。坐标为零起始 (行,位)，位 0 为 LSB。差分路径 {html.escape(b['differential_id'])}，移位 {b['differential_rotation']}；线性路径 {html.escape(b['linear_id'])}。</p><div class="scroll"><table><thead><tr><th>量</th><th>结果</th><th>含义</th></tr></thead><tbody>
<tr><td>原单条拼接 S</td><td>{number(b['score'])}</td><td>一条差分路径、中间连接及一对相同线性路径的评分</td></tr>
<tr><td>完整普通 DATF ρ̂</td><td>{number(datf)}</td><td>完整端点的普通 DATF 估计</td></tr>
<tr><td>HNW26 局部联合版 ρ̂</td><td>{number(full['correlation'])}</td><td>每列 1,024 个联合状态／差分相关分量</td></tr>
<tr><td>HNW26 零状态值掩码版 ρ̂</td><td>{number(zero['correlation'])}</td><td>每列 32 个差分相关分量</td></tr>
<tr><td>已有置换采样 ρ̂</td><td>{float(sample['rho']):.10f}</td><td>{int(sample['samples']):,} 对输入，95% 区间 [{lo:.10f}, {hi:.10f}]</td></tr></tbody></table></div><div class="callout"><p>{conclusion}</p></div><div class="states">{states('输入差分 Δ',b['input_difference'])}{states('输出掩码 Λ',b['output_mask'])}</div></section>
<section><h2>首轮分支与零值来源</h2><p>首个 χ 的全部 {full['completed_branches']} 个可行 DDT 分支均已计算，概率和严格为 1。{cause}</p><div class="scroll"><table><thead><tr><th>分支</th><th>首轮输出差分位置</th><th>概率</th><th>联合版条件 ρ̂</th><th>零值版条件 ρ̂</th><th>联合版完全退化轮次</th><th>加权贡献</th></tr></thead><tbody>{''.join(branch_rows)}</tbody></table></div><p>“完全退化”指全部局部相关向量只保留零掩码分量 1，其余分量严格为 0；轮次从所评估区间的首个 χ 开始计数。使用 GMP 任意精度整数和 2 的幂分母，不使用浮点零阈值。</p></section>
<section><h2>计算设置与解释范围</h2><p>首轮精确枚举 DDT 分支；之后逐列传播局部相关向量。联合版保留条件状态／差分联合分布与后续轮常数，简化版只保留零状态值掩码切片。两列输出在每个首轮分支内采用局部相关分量乘积读出，再按分支概率求和。这一多列读出同样依赖列间独立近似。</p><p>HNW 的核心近似是将跨列联合相关度替换为局部量的乘积。即使保留每列 1,024 个联合分量，也没有保留不同列之间的全部联合关系。本次结果只评价这个端点和上述 HNW 设置，不说明所有 HNW 分段或高阶块近似都必然得到相同结果。</p><p>输出是完整端点的模型相关度，没有再乘原 1+1+2 路径的 p<sub>D</sub> 或 c<sub>L</sub>²。HNW 对应偏差 ε̂ = ρ̂/2，联合版为 {number(full['bias'])}，简化版为 {number(zero['bias'])}。已有采样直接沿用，未重新采样。</p></section>
<section><h2>验证与文件</h2><p>算术后端自测通过；逐分支重新核对 DDT 概率、概率和、条件相关度与加权贡献和；另以独立 Python Fraction 递归复算简化版，得到 {html.escape(str(oracle))}。</p><p><a href="results.json">汇总 JSON</a> · <a href="configuration.json">端点与设置</a> · <a href="validation.json">独立核验</a> · <a href="backend_validation.json">后端自测</a> · <a href="jobs/{eid}.full.json">联合版全部分支</a> · <a href="jobs/{eid}.zero.json">简化版全部分支</a> · <a href="../linear_two_matches.html">返回匹配报告</a> · <a href="../sampling.csv">已有采样计数</a></p></section></main></body></html>'''
    (out/'hnw26_endpoint.html').write_text(page)
    save(out/'reproducibility.json',{'files':[{'path':str(f.resolve()),'sha256':hashlib.sha256(f.read_bytes()).hexdigest()}
        for f in sources+[f for f in out.rglob('*') if f.is_file() and f.name!='reproducibility.json']]})
    print('REPORT',out/'hnw26_endpoint.html',flush=True)


if __name__=='__main__':main()
