#!/usr/bin/env python3
"""Render the strict 1+2+1 -> 1+2+2 extension and its ordinary DATF result."""
import gzip
import hashlib
import html
import json
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlsplit

from extend_121_tail import ROOT, decode
from hnw26_compare import save


def state(title, words):
    text = '\n'.join(f'x{i}: {int(w,16):016x}' for i,w in enumerate(words))
    return f'<div><h3>{html.escape(title)}</h3><pre>{text}</pre></div>'


def main():
    root = ROOT / 'results/extend_121_tail'
    data = json.loads((root/'results.json').read_text())
    config = json.loads((root/'configuration.json').read_text())
    validation = json.loads((root/'validation.json').read_text())
    for path, digest in config['source_hashes'].items():
        assert hashlib.sha256(Path(path).read_bytes()).hexdigest() == digest
    strict = next(r for r in data['results'] if r['id'] == 'strict_append')
    best = next(r for r in data['results'] if r['id'] == 'reoptimized_tail')
    e = strict['extension']
    with gzip.open(strict['DATF']['artifact'], 'rt') as f:
        datf = json.load(f)['result']
    assert datf['status'] == 'completed' and datf['factor_audit']['status'] == 'PASS'
    assert datf['correlation']['numerator'] == '0'
    assert e['added_round_active_columns'] == len(e['second_input']['active_columns']) == 62
    assert e['first_output']['words'] == ['0x'+w for w in data['source_candidate']['output_mask']]
    trace_rows = ''.join(f'<tr><td>{f["step"]}</td><td>{f["variables"]}</td><td>{f["terms"]}</td>'
                         f'<td>{"2<sup>−128</sup>" if f["step"] == 0 else "0"}</td>'
                         f'<td>{"独立二次型精确求和" if f["step"] == 0 else "均匀平移翻转证书"}</td></tr>'
                         for f in datf['factor_audit']['factors'])
    all_records = json.loads((root/'all_extensions.json').read_text())['records']
    scan_rows = ''.join(f'<tr><td>{r["rank"]}</td><td>{"原掩码" if r["preserves_original_last_mask"] else "调整输出"}</td>'
                        f'<td>{r["first_round_squared_weight"]}</td><td>{r["added_round_active_columns"]}</td>'
                        f'<td>{r["added_round_squared_weight"]}</td><td>{r["tail_squared_weight"]}</td></tr>' for r in all_records)
    page = f'''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Gaston 1+2+1 向后延长一轮 · 普通 DATF</title><style>
body{{margin:0;background:#f1f5f7;color:#17354b;font:16px/1.85 system-ui,sans-serif}}main{{max-width:1120px;margin:auto;padding:34px 24px 72px}}h1{{font-size:32px;line-height:1.4}}h2{{font-size:23px}}h3{{font-size:16px}}section{{background:white;border:1px solid #d5e2e9;border-radius:12px;padding:25px;margin:24px 0}}a{{color:#146c96}}.cards{{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:16px}}.card{{padding:18px 22px;background:white;border:1px solid #d5e2e9;border-top:4px solid #307c9b;border-radius:8px}}.card strong{{display:block;font-size:35px;line-height:1.5}}.card span{{font-size:14px;color:#56707f}}.states{{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:20px}}pre{{background:#f0f5f7;padding:15px;overflow:auto;font:14px/1.85 ui-monospace,monospace}}.callout{{background:#fff4df;border-left:4px solid #ac761e;padding:16px 20px}}.note{{background:#edf5f8;padding:16px 20px}}table{{border-collapse:collapse;width:100%;font-size:14px}}td,th{{padding:12px 10px;border-bottom:1px solid #d5e2e9;text-align:left}}th{{background:#edf4f7}}.scroll{{overflow:auto}}.formula{{font:19px/1.8 Georgia,serif;text-align:center;overflow:auto;padding:16px;background:#f0f5f7}}.muted{{font-size:14px;color:#536c7b}}sup{{line-height:0}}details{{margin:20px 0}}summary{{cursor:pointer;color:#146c96}}@media(max-width:740px){{main{{padding:20px 14px 55px}}h1{{font-size:27px}}section{{padding:18px}}.states,.cards{{grid-template-columns:1fr}}}}
</style></head><body><main><p class="muted">GASTON · 1+2+1 → 1+2+2 · 2026-10-17</p><h1>固定原路径，线性部分向后延长一轮</h1>
<p>固定原 1+2+1 候选的输入差分、差分段、中间段、原线性轮输入及输出掩码，在末端补上 L、轮常数和一层 χ。选取新增 χ 层局部相关度最大的一个输出掩码，形成五轮 1+2+2 候选，并计算完整普通 DATF。</p>
<div class="cards"><div class="card"><span>新增 χ 活跃列数</span><strong>62</strong><span>已按逐列掩码精确核对</span></div><div class="card"><span>新增 χ 最优 c²</span><strong>2<sup>−128</sup></strong><span>对应 |c| = 2<sup>−64</sup></span></div><div class="card"><span>完整五轮普通 DATF</span><strong>0</strong><span>计算完成，独立因子复核通过</span></div></div>
<p class="callout">此前对话中将新增活跃列数说成 61，准确计数为 62。新增 χ 输入掩码有 170 个活跃比特，分布于 62 列；第 30、47 列不活跃。局部权重与 DATF 结果不受这处口头计数修正影响。</p>
<section><h2>延长得到的路径</h2><p>输入差分 Δ = &#123;(0,1)&#125;。原线性轮输入与输出掩码均为 &#123;(1,0),(3,54)&#125;。坐标使用从 0 起的 (行,位)，位 0 是 LSB。</p>
<div class="formula">a₁ → χ → b₁ → L → a₂ → χ → b₂，&nbsp; b₁ = Lᵀa₂，&nbsp; a₂ = L<sup>−T</sup>b₁</div>
<p>掩码向后延长需要求 L 的逆转置。固定 b₁ 后，a₂ 唯一确定；这里没有将正向状态传播 L(b₁) 当作掩码传播。计算包含 ρ<sub>east</sub>、θ、ρ<sub>west</sub> 的完整作用。</p>
<div class="states">{state('原线性轮输出 b₁，保持不变', e['first_output']['words'])}{state('新增 χ 的输入掩码 a₂', e['second_input']['words'])}{state('最终输出掩码 b₂，本次 DATF 使用', e['final_output']['words'])}{state('完整五轮输入差分 Δ', data['source_candidate']['input_difference'])}</div>
<div class="scroll"><table><thead><tr><th>量</th><th>数值</th><th>范围</th></tr></thead><tbody>
<tr><td>原差分路径概率 p<sub>D</sub></td><td>2<sup>−2</sup></td><td>首轮 χ 的指定 DDT 转移</td></tr>
<tr><td>两轮中间段普通 DATF</td><td>+2<sup>−2</sup></td><td>已按五轮区间的常数复算</td></tr>
<tr><td>原线性轮 c₁²</td><td>2<sup>−4</sup></td><td>保持原线性转移</td></tr>
<tr><td>新增线性轮 c₂²</td><td>2<sup>−128</sup></td><td>固定 a₂ 下的最优局部 LAT 乘积</td></tr>
<tr><td>两轮线性尾段 c<sub>L</sub>²</td><td>2<sup>−132</sup></td><td>单条线性路径</td></tr>
<tr><td>1+2+2 拼接评分 S</td><td>+2<sup>−136</sup></td><td>p<sub>D</sub> × ρ̂<sub>middle</sub> × c<sub>L</sub>²</td></tr>
<tr><td>完整五轮普通 DATF ρ̂</td><td><strong>0</strong></td><td>直接评估 Δ 与 b₂ 的完整端点</td></tr>
</tbody></table></div><p>本次所选两轮线性路径的带符号相关度为 −2<sup>−66</sup>，其平方为 2<sup>−132</sup>。共有 2<sup>26</sup> 个同样达到新增 χ 最优局部权重的末端掩码；本次选取每列最小比特重量、再按五比特整数值取最小的一条。完整 DATF 只计算了这个代表，没有穷尽所有同分末端掩码。</p></section>
<section><h2>普通 DATF 如何得到 0</h2><p>直接构造五轮输出差分的掩码奇偶多项式，初始包含 {datf['target_terms']} 项、{datf['target_variables']} 个 DAG 变量。采用普通 DATF Algorithm 1，以精确二次型算法计算 Bias；未使用增强 DATF。精确二次型计算只是加速原有 Bias 子过程，不改变普通 DATF 的独立性假设。</p>
<div class="scroll"><table><thead><tr><th>顶层消元步骤</th><th>因子变量数</th><th>因子项数</th><th>因子相关度</th><th>独立复核</th></tr></thead><tbody>{trace_rows}</tbody></table></div>
<p>第一次提取的因子为 2<sup>−128</sup>；继续替换孤立变量后，下一因子在模型的独立边缘分布下相关度为 0，故整体模型乘积为 0。零因子存在平移证书 f(z⊕h)=f(z)⊕1，正负值在该模型分布下配对抵消。</p>
<p>该证书核对的是 DATF 因子的算术，不证明真实置换相关度为零。32 组完整 DAG 与直接置换求值核对全部通过，因子独立复核为 PASS；计算正常完成，没有把超时、资源限制或浮点小值记成 0。</p>
<p class="note">拼接评分 +2<sup>−136</sup> 与完整 DATF 的 0 是不同层次的量。前者只使用指定路径和中间模型值，后者从最终完整端点重新执行普通 DATF；二者都不能直接当作真实五轮相关度。</p></section>
<section><h2>轮常数与完整区间</h2><p>沿用此前“取参考实现末 R 轮”的约定，五轮对应参考第 8–12 轮；旧四轮掩码结构位于第 8–11 轮，新 χ 位于第 12 轮。仅省略整个区间开头的第一个 L，R=5 含 5 个 χ 和 4 个 L。</p>
<p>首轮常数吸收到均匀输入，后续常数依次为 0x78、0x69、0x5a、0x4b。两轮中间段的内部常数由原来的 0x5a 变为 0x69；两种设置均复算得到 +1/4。另核对了平移后的旧四轮完整端点，普通 DATF 仍为 3/2048。</p>
<p>补充计算的完整 HNW26 联合版与零状态值掩码版均为 0。它们逐一覆盖首轮全部 4 个 DDT 分支；这些结果同样仅是局部列近似结论。</p></section>
<details><summary>补充：保持线性段入口不变，重选两轮尾段的全部 100 个首轮输出</summary><section><h2>固定入口的两轮线性路径最优化</h2><p>两个活跃列各有 10 个非零 LAT 输出，共 100 种原线性轮输出。每一种都解出唯一 a₂，再逐列求新增 χ 的最大绝对 LAT 值，因此单路径权重优化覆盖全部最终掩码，无需逐个枚举其笛卡尔积。</p>
<p>最好的尾段 c<sub>L</sub>² 为 2<sup>−124</sup>，两轮的平方权重分别为 6、118，新增 χ 有 58 个活跃列。对应拼接评分 +2<sup>−128</sup>。这项改变了原线性轮输出，属于固定入口的重新优化；它没有替代上方的严格延长结果。对其所选代表，完整普通 DATF 和两种 HNW26 均为 0。</p>
<div class="scroll"><table><thead><tr><th>排名</th><th>是否保持原输出</th><th>第一轮平方权重</th><th>新增活跃列数</th><th>新增平方权重</th><th>尾段平方权重</th></tr></thead><tbody>{scan_rows}</tbody></table></div>
<p>这证明的是固定线性入口下的最佳单路径权重，未证明完整线性包或完整差分线性相关度的最优性。</p></section></details>
<section><h2>结果与复现文件</h2><p><a href="results.json">结果 JSON</a> · <a href="configuration.json">固定端点与配置</a> · <a href="validation.json">验证记录</a> · <a href="datf_jobs/strict_append.json.gz">严格延长的完整 DATF 多项式、DAG 与轨迹</a> · <a href="all_extensions.json">全部 100 种尾段的最优权重与掩码</a> · <a href="../datf_pools/search_final/datf_pool_matches.html">原 1+2+1 报告</a></p>
<pre>python3 tools/extend_121_tail.py
python3 tools/extend_121_report.py</pre>
<p>线性层逆转置矩阵秩为 320；全部 100 个桥接均用正向 L 的 320 个基向量独立核对，共 32,000 项。DATF 中间段另由 C++ 与 Python 两套实现交叉核对。</p></section></main></body></html>'''
    path = root/'extend_121_tail.html'
    path.write_text(page)
    class Links(HTMLParser):
        def __init__(self):
            super().__init__()
            self.links = []
        def handle_starttag(self, tag, attrs):
            attrs = dict(attrs)
            if 'href' in attrs:
                self.links.append(attrs['href'])
    parsed = Links()
    parsed.feed(page)
    for link in parsed.links:
        u = urlsplit(link)
        if not u.scheme and u.path:
            assert (root/unquote(u.path)).exists(), link
    save(root/'report_validation.json', {'status':'PASS', 'local_links_checked':len(parsed.links),
         'strict_added_active_columns':62, 'inactive_columns':[j for j in range(64) if j not in e['second_input']['active_columns']],
         'source_hashes_verified':True, 'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'html_sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
    print(path)


if __name__ == '__main__':
    main()
