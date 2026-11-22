#!/usr/bin/env python3
"""Render an offline HNW26 / DATF / sampling comparison report."""
import argparse
import hashlib
import html
import json
import os
from pathlib import Path

from datf_merged_report import CSS
from datf_pool_report import number,state_block
from hnw26_compare import ROOT,fraction,save


def main():
    p=argparse.ArgumentParser();p.add_argument('--results',type=Path,default=ROOT/'results/hnw26_comparison')
    root=p.parse_args().results.resolve()
    data=json.loads((root/'results.json').read_text());validation=json.loads((root/'validation.json').read_text())
    assert data['status']=='completed' and validation['status']=='PASS'
    rows=data['results'];lookup={r['id']:r for r in rows}
    def esc(x):return html.escape(str(x),quote=True)
    def link(path,label):return '<a href="'+esc(os.path.relpath(path,root))+'">'+esc(label)+'</a>'
    def val(r,mode):return number(r['hnw'][mode]['correlation'])
    def histogram(r,mode='full'):
        return '；'.join(f'第 {k} 轮：{v:,} 支' for k,v in r['hnw'][mode]['absorbed_histogram'].items()) or '未完全退化'
    selected_rows=[];table=[];details=[]
    for i,r in enumerate(rows):
        b=r['candidate'];full=r['hnw']['full'];zero=r['hnw']['zero']
        agree=fraction(r['datf_correlation'])==fraction(full['correlation'])
        category='原 DATF 非零' if r['category']=='DATF_nonzero' else '原 DATF 为零'
        cells=[f'<b>{b["R"]} 轮</b><small>来源 {r["original_split"]}<br>HNW {r["hnw_split"]}</small>',
               f'<span class="tail">{esc(r["id"])}</span>',number(b['score']),number(r['datf_correlation']),
               val(r,'full'),val(r,'zero'),f'{full["completed_branches"]:,}',
               f'<span class="badge {"good" if agree else "pending"}">{"一致" if agree else "不同"}</span>',
               f'<a href="#case-{i}">边界与分支</a>']
        table.append(f'<tr data-category="{r["category"]}" data-agree="{str(agree).lower()}">'+''.join('<td>'+x+'</td>' for x in cells)+'</tr>')
        if r['category']=='DATF_nonzero':
            selected_rows.append('<tr>'+''.join('<td>'+x+'</td>' for x in [r['original_split'],str(b['R']),number(r['datf_correlation']),val(r,'full'),val(r,'zero')])+'</tr>')
        state=''.join(state_block(label,b[key]) for label,key in [('完整输入差分 Δ','input_difference'),('完整输出掩码 Λ','output_mask')])
        old_formula='pD = '+number(b['differential_probability'])+'；ρ̂M = '+number(b['middle_DATF_correlation'])+'；cL² = '+number(b['linear_squared_correlation'])+'；S = '+number(b['score'])
        details.append(f'<details id="case-{i}"><summary>{r["original_split"]} · {esc(r["id"])} · {category}</summary><p>选取该类原拼接排名第 1 的端点。差分路径 <code>{esc(b["differential_id"])}</code>，线性路径 <code>{esc(b["linear_id"])}</code>，差分共同左移 {b["differential_rotation"]} 位；采用参考实现末 {b["R"]} 轮，零起始轮号 {b["start_round_index"]}。</p><div class="states">{state}</div><p>原路径分段：{old_formula}。这些原路径权重仅供对照，没有再次乘入 HNW 结果。</p><p>HNW 重评采用 {r["hnw_split"]}，首轮完整枚举 {full["completed_branches"]:,} 个 DDT 输出差分分支，概率和严格为 1。每个分支内传播完整区间并计算最终掩码，最后进行带符号加权。</p><div class="formula">DATF ρ̂ = {number(r["datf_correlation"])}；HNW 局部联合版 ρ̂ = {val(r,"full")}，ε̂ = {number(full["bias"])}；HNW 零值掩码版 ρ̂ = {val(r,"zero")}</div><p>局部向量完全退化为仅零掩码分量为 1 的分支：{histogram(r)}。在目标轮数处有 {full["nonzero_branch_correlations"]} 个分支的最终模型相关度非零。</p><p>{link(Path(full["artifact"]),"局部联合版：全部分支与精确贡献")} · {link(Path(zero["artifact"]),"零值掩码版：全部分支与精确贡献")}</p></details>')
    sampling=[]
    for e in data['existing_experiments']:
        if e['id'] not in lookup:continue
        r=lookup[e['id']];value=float(fraction(r['hnw']['full']['correlation']))
        inside=float(e['ci95_low'])<=value<=float(e['ci95_high'])
        sampling.append('<tr>'+''.join('<td>'+x+'</td>' for x in [esc(e['id']),
            '独立确认' if e['source_file']=='confirmation.csv' else '首批',f'{int(e["samples"]):,}',
            f'{float(e["rho"]):.9f}',f'[{float(e["ci95_low"]):.9f}, {float(e["ci95_high"]):.9f}]',
            number(r['datf_correlation']),val(r,'full'),'是' if inside else '否'])+'</tr>')
    thunder=Path('/home/yxb/ngcc/spec/hash-33-spec_thunder.pdf')
    page=f'''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Gaston · HNW26 与 DATF 固定端点对照</title><style>{CSS}
#comparison{{min-width:1200px}}.equation{{font-family:ui-monospace,monospace;line-height:1.9;padding:14px;background:#eef3f0;border-radius:6px;overflow:auto}}.proof{{border-left:3px solid #779385;padding-left:18px}}.highlight{{color:#925b15}}.compact td{{padding:9px 12px}}</style></head><body><main>
<header><div class="eyebrow">GASTON / HNW26 × ORDINARY DATF / 2026-10-16</div><h1>固定端点的 HNW26 逐轮近似复评</h1><p class="subtitle">保留已有输入差分与最终输出掩码，使用 Thunder 的 r<sub>d</sub>=1、r<sub>l</sub>=0 设置重评完整区间。并列计算 HNW 的局部联合信息版与零状态值掩码简化版，使用精确整数算术；仍省略整个区间的第一个线性层。</p><nav class="nav"><a href="#findings">三个非零候选</a><a href="#results">全部结果</a><a href="#sampling">已有采样</a><a href="#method">方法与范围</a><a href="#why-zero">零值原因</a><a href="#details">端点边界</a><a href="#files">验证与文件</a></nav></header>
<div class="cards"><div class="card"><span>固定端点</span><strong>{len(rows)}</strong><span>3 个 DATF 非零 + {len(rows)-3} 个零例代表</span></div><div class="card"><span>HNW 版本</span><strong>2 版</strong><span>每列 1,024 / 32 个相关分量</span></div><div class="card"><span>HNW 非零端点</span><strong>1</strong><span>1+1+1 的三轮候选</span></div><div class="card"><span>计算用时</span><strong>{data['elapsed_seconds']:.1f} 秒</strong><span>{data['workers']} 个原生进程槽位；不含核验</span></div></div>
<section id="findings"><h2>三个已有 DATF 非零端点的结果</h2><div class="scroll"><table class="compact"><thead><tr><th>原搜索分段</th><th>总轮数</th><th>完整普通 DATF</th><th>HNW 局部联合版</th><th>HNW 零值掩码版</th></tr></thead><tbody>{''.join(selected_rows)}</tbody></table></div><div class="callout"><b>本次两版 HNW 的结果逐项相同：三轮为 −2<sup>−4</sup>，其余选定端点均为 0。</b><p>两个四轮非零 DATF 端点在 HNW 下都退化为模型零。其中 1+2+1 已有采样信号，因此这里的 HNW 零不能解释为真实相关度为零；对于尚无采样的 1+1+2，模型分歧本身不能决定哪一个估计更接近真实值。</p></div><p>这是一组固定端点对照，没有重新寻找输入差分或输出掩码。结论只覆盖本报告所列端点和当前逐列近似设置。</p></section>
<section id="results"><h2>26 个固定端点的对照总表</h2><p>“每一类”按综合报告每行定义，即差分池 × 中间轮数 × 线性尾段。对 22 行论文尾段的全零类及单轮线性池的 2+1+1 零例，各取原拼接首名，共 23 个零例。另加入指定的 3 个非零例。</p><div class="filters"><label>原 DATF <select id="category"><option value="">全部</option><option value="DATF_nonzero">非零候选</option><option value="DATF_zero">零例代表</option></select></label><label>与 HNW 比较 <select id="agree"><option value="">全部</option><option value="false">结果不同</option><option value="true">结果一致</option></select></label><span id="shown" aria-live="polite"></span></div><div class="scroll"><table id="comparison"><thead><tr><th>轮数 / 分段</th><th>候选类别</th><th>原拼接 S</th><th>完整 DATF ρ̂</th><th>HNW 局部联合 ρ̂</th><th>HNW 零值掩码 ρ̂</th><th>首轮分支数</th><th>模型对照</th><th>详情</th></tr></thead><tbody>{''.join(table)}</tbody></table></div><p class="legend">原分段是端点的搜索来源；HNW 分段统一为 1+(R−1)+0。例如原 2+1+2 五轮候选，在 HNW 中按 1+4+0 计算。HNW 结果为完整端点的模型相关度，不是把原路径概率与新中间段估计再次相乘。</p></section>
<section id="sampling"><h2>与已有完整置换采样对照</h2><p>保留已有实验数据，没有新增采样。区间为各项的 95% Wilson 区间。</p><div class="scroll"><table id="samples"><thead><tr><th>端点</th><th>实验</th><th>样本对数</th><th>经验相关度</th><th>95% 区间</th><th>DATF</th><th>HNW</th><th>HNW 在区间内</th></tr></thead><tbody>{''.join(sampling)}</tbody></table></div><p><b>1+2+1 四轮候选</b>的独立确认采样为 0.004008472，区间 [0.003769220, 0.004247723]，排除 0。HNW 给出 0；DATF 给出 3/2048 ≈ 0.001464844，也在区间外。因此 DATF 保留了非零信号，但两种模型都没有准确预测该端点的量级。</p><p class="legend">三轮候选的模型值 −1/16 与现有采样相容；2+1+1 的区间包含 0，只表示目前样本未分辨出非零信号，不构成零相关度证明。模拟使用原参考轮常数和首个 L 省略的相同边界。</p></section>
<section id="method"><h2>使用的 HNW26 方法与实现选择</h2><p>{link(thunder,'Thunder 规范')}第 4.3.7 节及 Table 13（第 27–28 页）明确使用首轮精确差分分布和 r<sub>d</sub>=1、r<sub>l</sub>=0。该节没有说明采用哪一种局部向量版本，也未提供对应实验程序；本次据 <a href="https://github.com/hukaisdu/Round-based-approximation/blob/783e8edc560b848d0f2ba63f11838aecd23da40f/Ascon/ascon1st.cpp">HNW 作者公开实现</a>独立实现两种更新，替换为 Gaston 的 χ 与线性层，并明确记录多活跃输出的读出方法。不是对 Thunder 私有实验程序的逐行复现。</p>
<h3>首轮精确枚举，后续逐列近似</h3><p>对完整输入差分 Δ 的所有活跃 χ 列，枚举首轮全部输出差分 β，概率为各列 DDT[Δj,βj]/32 的乘积。局部联合版保留首轮分支条件下 χ(X) 的实际局部分布；零值掩码版仅保留固定差分。每个分支独立传播到最后，再按首轮概率对带符号相关度求和。</p>
<div class="equation">ρ̂(Δ,Λ) = Σ<sub>β</sub> Pr[χ(X) ⊕ χ(X⊕Δ)=β] · ρ̂<sub>β</sub>(Λ)</div><p>即使端点来源于两轮 6→6 路径，本次也没有固定那条中间路径。六个单比特活跃列产生 4⁶=4,096 个首轮分支，均完整计算。</p>
<h3>每列局部联合信息版：1,024 维</h3><div class="equation">q<sub>j</sub>(a,b) = E[ (−1)<sup>a·Xj ⊕ b·ΔXj</sup> ]<br>q′<sub>j</sub>(a,b) = Σ<sub>u,v</sub> C[a⊕b,u⊕v] C[b,v] q<sub>j</sub>(u,v)<br>q<sub>L,j</sub>(a,b) ≈ ∏<sub>k</sub> q<sub>k</sub>((Lᵀ embed(j,a))<sub>k</sub>, (Lᵀ embed(j,b))<sub>k</sub>)</div><p>C 是归一化的 χ 线性相关矩阵。非线性更新用 Walsh 变换和联合状态置换精确实现，并与独立的双 LAT 公式逐项交叉核对。线性更新中的跨列乘积才是这一层的近似。</p>
<h3>零状态值掩码简化版：32 维</h3><div class="equation">q<sub>j</sub>(b) = q<sub>j</sub>(0,b)<br>q′<sub>j</sub>(b) = Σ<sub>v</sub> C[b,v]² q<sub>j</sub>(v)<br>q<sub>L,j</sub>(b) ≈ ∏<sub>k</sub> q<sub>k</sub>((Lᵀ embed(j,b))<sub>k</sub>)</div><p>多列输出掩码在每个首轮分支内部读作 ∏<sub>j</sub>q<sub>j</sub>(0,Λj)，然后才进行分支加权；没有先对各列混合分支再取乘积。该读出也依赖局部列分解，超出 Thunder Table 13 的单活跃输出枚举范围。</p><p>完整区间仍为 χ ∘ L ∘ … ∘ L ∘ χ，采用参考实现末 R 轮。局部联合版保留后续参考轮常数，零值掩码版在模型中不受常数影响；下面的结构核验进一步证明，在本批端点的局部联合模型中，常数也不会改变结果。因此保留与省略后续常数在这批模型计算上等价，不能推广为真实置换的常数无关性。</p></section>
<section id="why-zero"><h2>为什么两版相同，而且四轮端点均变为零</h2><div class="proof"><p><b>第一步：首个 L 消去局部非零状态值掩码分量。</b>对本批每个不同输入差分、每个非零局部状态值掩码，都核对了 Lᵀ 的支撑必然触及首轮输入的某个不活跃列。该列的条件基值仍均匀，给乘积近似提供一个零因子。</p><p><b>第二步：后续 L 重复消去这些分量。</b>穷举核对全部 64×31×32=63,488 个局部掩码对 (a≠0,b)：Lᵀ 的某个输入列总满足 u₀≠0，且 u₁=0 或 u₁=u₀。对应项分别是原状态或配对状态的非平凡线性掩码；当两者边缘均匀时，该项为 0。χ 的双射性保持两边状态的均匀边缘。</p><p>所以每次进入 χ 之前，局部联合版只剩 a=0 的分量，其更新化为零值掩码版的 LAT 平方更新。两版相同有本模型内的结构依据，不只是数值碰巧相等。</p></div><p>对两个 DATF 非零的四轮端点，4 个首轮分支中有 2 支在第 3 轮退化，另 2 支在第 4 轮退化。这里“退化”指近似局部向量仅保留零掩码分量 1，所有其他分量都严格为 0。全部 25 个选定的四轮及以上端点，每一支最终都为 0；此次没有发生正负非零分支抵消。</p><div class="callout">这些零由任意精度整数与 2 的幂分母精确计算，排除了浮点下溢。它们仍是局部独立性模型中的零；1+2+1 的现有采样已经表明，这种逐列近似在该端点丢失了真实的跨列相关信息。不能把本结论写成 Gaston 四轮真实相关度为零。</div><p>{link(root/'model_equivalence.json','模型内等价性核验与推导')} · {link(root/'linear_mask_witnesses.json.gz','全部 63,488 个掩码对的零因子见证')}</p></section>
<section id="details"><h2>逐项端点与精确分支记录</h2>{''.join(details)}</section>
<section id="files"><h2>计算验证与复现</h2><p>后端使用 GMP 任意精度整数，不使用浮点零阈值。验证包含 4,096 个双 LAT / Walsh 更新坐标、317 个可行条件 DDT 分支、1,024 个零掩码切片坐标、65,536 个线性传播坐标及同样数量的常数相位坐标；另与旧 Table 13 实现核对，并用独立 Python Fraction 递归计算三个指定端点。</p><p>对本次 52 次原生计算，逐项核对首轮分支去重、DDT 精确概率、概率和、带符号贡献和及端点参数。结果没有未知、截断或资源上限代替的零。</p><pre>make build/hnw_endpoint
./build/hnw_endpoint --selftest &gt; results/hnw26_comparison/backend_validation.json
python3 tools/hnw26_compare.py --workers 24
python3 tools/hnw26_audit.py
python3 tools/hnw26_report.py</pre><p class="legend">重跑默认复用配置完全相同的已完成端点；更换设置时使用新 --output 目录。HNW 主程序额外需要本机已安装的 GMP 开发库（-lgmpxx -lgmp）。</p><div class="links">{link(root/'results.json','完整结果汇总')}{link(root/'configuration.json','端点清单与设置')}{link(root/'validation.json','独立核验记录')}{link(root/'backend_validation.json','算术后端核验')}{link(root/'reproducibility.json','代码与结果 SHA-256')}{link(ROOT/'src/hnw_endpoint.cpp','C++ 精确计算器')}{link(ROOT/'results/datf_pools/merged/datf_merged.html','原综合报告')}</div><p>方法来源：<a href="https://eprint.iacr.org/2026/358">Hu–Niu–Wang，Round-Based Approximation of (Higher-Order) Differential-Linear Correlation</a>；作者公开代码固定提交 <code>783e8edc560b848d0f2ba63f11838aecd23da40f</code>。Thunder 参数来源见第 4.3.7 节。</p></section><footer>研究对象：Gaston 减轮置换 · 原始普通 DATF 与 HNW26 局部逐轮近似 · 首个 L 省略 · 2026-10-16</footer></main>
<script>function update(){{let c=document.querySelector('#category').value,a=document.querySelector('#agree').value,n=0;document.querySelectorAll('#comparison tbody tr').forEach(r=>{{r.hidden=!!((c&&r.dataset.category!==c)||(a&&r.dataset.agree!==a));if(!r.hidden)n++}});document.querySelector('#shown').textContent=`显示 ${{n}} / 26 行`}}document.querySelectorAll('select').forEach(e=>e.addEventListener('change',update));update();function openHash(){{let el=document.getElementById(location.hash.slice(1));if(el&&el.tagName==='DETAILS')el.open=true}}window.addEventListener('hashchange',openHash);openHash();document.querySelectorAll('a[href^="#case-"]').forEach(a=>a.addEventListener('click',()=>{{document.getElementById(a.hash.slice(1)).open=true}}));</script></body></html>'''
    (root/'hnw26_comparison.html').write_text(page)
    source=[ROOT/name for name in ['src/hnw_endpoint.cpp','include/gaston.hpp','include/estimator.hpp','tools/hnw26_compare.py','tools/hnw26_audit.py','tools/hnw26_report.py','Makefile']]
    source+=[thunder,ROOT/'results/datf_pools/merged/merged_results.json']
    files=source+[f for f in root.rglob('*') if f.is_file() and f.name!='reproducibility.json']
    save(root/'reproducibility.json',{'scope':'HNW26 fixed-endpoint comparison: source code, parameters, exact branch results and report',
        'files':[{'path':str(f.resolve()),'bytes':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()} for f in files]})
    print(root/'hnw26_comparison.html')


if __name__=='__main__':main()
