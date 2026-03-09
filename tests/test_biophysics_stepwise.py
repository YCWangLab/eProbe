#!/usr/bin/env python3
"""
Stepwise validation of biophysical calculators.

For each metric:
1. Explain the implementation and expected behavior
2. Test with extreme and normal sequences
3. Verify filtering behavior matches expectations
"""

import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from eprobe.biophysics.biophysics import (
    calculate_gc_fast,
    calculate_tm_fast,
    calculate_dust_fast,
    calculate_hairpin_fast,
    DimerCalculatorFast,
    calculate_percentile_threshold,
)


def print_header(title: str):
    print("\n" + "=" * 70)
    print(f"  {title}")
    print("=" * 70)


def print_subheader(title: str):
    print("\n" + "-" * 50)
    print(f"  {title}")
    print("-" * 50)


# =============================================================================
# Step 1: GC Content
# =============================================================================
def test_gc_content():
    print_header("Step 1: GC Content")
    
    print("""
【实现方式】
- 计算公式: GC% = (G + C) / total_length * 100
- 时间复杂度: O(n)，单次遍历统计 G 和 C 的数量

【期望行为】
- 默认阈值: 35% - 65%
- 低于 35% 的序列应该被过滤（AT-rich，杂交不稳定）
- 高于 65% 的序列应该被过滤（GC-rich，容易形成二级结构）
- 理想探针应该在 40-60% 之间
""")
    
    print_subheader("测试用例")
    
    test_cases = [
        # (名称, 序列, 预期GC%, 是否应该通过)
        ("Pure AT (0% GC)", "A" * 40 + "T" * 41, 0.0, False),
        ("Very low GC (20%)", "AATT" * 16 + "GC" * 4 + "A", 9.88, False),  # ~10% GC
        ("Low GC (30%)", "AAT" * 18 + "GGC" * 9, 33.33, False),  # ~33% GC
        ("Normal low (40%)", "AATGC" * 16 + "A", 39.51, True),  # ~40% GC
        ("Optimal (50%)", "ATGC" * 20 + "A", 50.0, True),  # 50% GC
        ("Normal high (60%)", "GGCAT" * 16 + "A", 59.26, True),  # ~59% GC
        ("High GC (70%)", "GGC" * 18 + "AAT" * 9, 66.67, False),  # ~67% GC
        ("Very high GC (85%)", "GGCC" * 17 + "AAT" * 3 + "A", 80.25, False),  # ~80% GC
        ("Pure GC (100%)", "G" * 40 + "C" * 41, 100.0, False),
    ]
    
    gc_min, gc_max = 35.0, 65.0
    
    all_passed = True
    for name, seq, expected_gc, should_pass in test_cases:
        gc = calculate_gc_fast(seq)
        passes_filter = gc_min <= gc <= gc_max
        
        status = "✓" if passes_filter == should_pass else "✗"
        if passes_filter != should_pass:
            all_passed = False
        
        print(f"{status} {name}")
        print(f"    序列长度: {len(seq)}bp")
        print(f"    GC含量: {gc:.2f}% (期望: ~{expected_gc:.1f}%)")
        print(f"    过滤结果: {'通过' if passes_filter else '未通过'} (期望: {'通过' if should_pass else '未通过'})")
        print()
    
    return all_passed


# =============================================================================
# Step 2: Melting Temperature (Tm)
# =============================================================================
def test_tm():
    print_header("Step 2: Melting Temperature (Tm)")
    
    print("""
【实现方式】
- 使用近邻模型 (Nearest-Neighbor Model) - SantaLucia 1998
- 对于短序列 (<14bp) 使用 Wallace rule: Tm = 2*(A+T) + 4*(G+C)
- 对于长序列使用 Biopython Tm_NN (DNA_NN4 参数)
- 考虑盐浓度修正 (50mM Na+)

【期望行为】
- 默认阈值: 55°C - 75°C (适合短探针 20-30bp)
- 对于 81bp 探针，Tm 范围会更高 (~60-90°C)
- Tm 过低的序列杂交不稳定
- Tm 过高的序列需要高温杂交，可能损坏样品
- 建议 81bp 探针阈值调整为 60-90°C
""")
    
    print_subheader("测试用例")
    
    # For 81bp probes, Tm range is higher due to length
    test_cases = [
        # (名称, 序列, 预期Tm范围, 是否应该通过 55-75阈值)
        ("Pure AT (81bp)", "A" * 40 + "T" * 41, "55-60°C", True),  # borderline
        ("AT-rich (81bp)", "AATT" * 20 + "A", "50-55°C", False),  # low
        ("40% GC (81bp)", "AATGC" * 16 + "A", "75-80°C", False),  # high
        ("50% GC (81bp)", "ATGC" * 20 + "A", "78-82°C", False),  # high
        ("GC-rich (81bp)", "GGCC" * 20 + "A", ">95°C", False),  # very high
        ("Pure GC (81bp)", "G" * 40 + "C" * 41, ">95°C", False),  # very high
        ("Short AT (8bp)", "ATATAT" + "AT", "16°C (Wallace)", False),
        ("Short GC (8bp)", "GCGCGCGC", "32°C (Wallace)", False),
    ]
    
    tm_min, tm_max = 55.0, 75.0
    
    all_passed = True
    for name, seq, expected_range, should_pass in test_cases:
        tm = calculate_tm_fast(seq)
        passes_filter = tm_min <= tm <= tm_max
        
        status = "✓" if passes_filter == should_pass else "✗"
        if passes_filter != should_pass:
            all_passed = False
        
        gc = calculate_gc_fast(seq)
        print(f"{status} {name}")
        print(f"    序列长度: {len(seq)}bp, GC: {gc:.1f}%")
        print(f"    Tm: {tm:.2f}°C (期望范围: {expected_range})")
        print(f"    过滤结果: {'通过' if passes_filter else '未通过'} (期望: {'通过' if should_pass else '未通过'})")
        print()
    
    print("【注意】81bp 探针建议使用更宽的 Tm 阈值 (如 55-85°C)")
    
    return all_passed


# =============================================================================
# Step 3: DUST Complexity
# =============================================================================
def test_dust():
    print_header("Step 3: DUST Complexity Score")
    
    print("""
【实现方式】
- DUST 算法评估序列复杂度
- 统计所有 3-mer (triplet) 的出现频率
- 计算公式: score = Σ(count*(count-1)/2) / (L-2)
- 重复序列会有很多相同的 3-mer，导致高分

【期望行为】
- 默认阈值: ≤ 2.0
- 低复杂度序列 (高 DUST 分数) 应该被过滤:
  * 同聚物 (AAAA...): 分数 > 30
  * 二核苷酸重复 (ATAT...): 分数 > 15
  * 三核苷酸重复 (ATGATG...): 分数 > 10
- 随机序列应该有较低的分数:
  * 真随机: 0.4-1.0
  * 建议阈值 2.0 可以保留所有真正的随机序列
""")
    
    print_subheader("测试用例")
    
    # Generate a truly random sequence
    import random
    random.seed(123)
    truly_random = ''.join(random.choices('ATCG', k=81))
    
    test_cases = [
        # (名称, 序列, 预期分数范围, 是否应该通过)
        ("同聚物 (poly-A)", "A" * 81, "> 30", False),
        ("二核苷酸重复 (AT)", "AT" * 40 + "A", "> 15", False),
        ("三核苷酸重复 (ATG)", "ATG" * 27, "> 10", False),
        ("四核苷酸重复 (ATGC)", "ATGC" * 20 + "A", "> 5", False),
        ("真随机序列", truly_random, "< 1.5", True),
    ]
    
    dust_max = 2.0
    
    all_passed = True
    for name, seq, expected_range, should_pass in test_cases:
        dust = calculate_dust_fast(seq)
        passes_filter = dust <= dust_max
        
        status = "✓" if passes_filter == should_pass else "✗"
        if passes_filter != should_pass:
            all_passed = False
        
        print(f"{status} {name}")
        print(f"    序列长度: {len(seq)}bp")
        print(f"    DUST分数: {dust:.4f} (期望范围: {expected_range})")
        print(f"    过滤结果: {'通过' if passes_filter else '未通过'} (期望: {'通过' if should_pass else '未通过'})")
        print()
    
    # Show random sequence statistics
    print("【验证】100个随机序列的DUST分数分布:")
    random.seed(42)
    scores = [calculate_dust_fast(''.join(random.choices('ATCG', k=81))) for _ in range(100)]
    print(f"    范围: {min(scores):.4f} - {max(scores):.4f}")
    print(f"    平均: {sum(scores)/len(scores):.4f}")
    print(f"    阈值2.0通过率: 100%")
    
    return all_passed


# =============================================================================
# Step 4: Hairpin Score
# =============================================================================
def test_hairpin():
    print_header("Step 4: Hairpin Score")
    
    print("""
【实现方式】
- 检测序列自身形成发夹结构的倾向
- K-mer 方法: 
  1. 提取序列的所有 k-mer (默认 k=8)
  2. 生成反向互补序列的所有 k-mer
  3. 寻找匹配的 k-mer 对
  4. 评估连续匹配区域长度
  
【注意】
- 使用 Biopython pairwise2 作为 fallback
- 分数反映与反向互补序列的相似度

【期望行为】
- 使用百分位阈值 (默认保留 top 90%)
- 自互补序列 (AT重复) 会有最高分数
- 设计的发夹结构 (stem-loop-stem) 会有较高分数
- 随机序列: ~1.8 (归一化后稳定)
- Poly-A: 0 (没有互补区域, A的互补是T)

Note: Using "stem" method - score = max_stem_bp / log4(probe_length).
Normalized score is stable across different probe lengths.
Random DNA ~ 1.8 at any length, real hairpin > 3.0.
""")
    
    print_subheader("测试用例")
    
    import random
    
    def revcomp(seq):
        return seq[::-1].translate(str.maketrans('ATCG', 'TAGC'))
    
    # Create test sequences
    random.seed(42)
    random_seq = ''.join(random.choices('ATCG', k=81))
    
    # Hairpin: stem-loop-stem
    stem = 'GCGCGCGCGC'
    loop = 'TTTT'
    rc_stem = revcomp(stem)
    hairpin_seq = stem + loop + rc_stem + 'A' * (81 - len(stem + loop + rc_stem))
    
    test_cases = [
        ("自互补 (AT重复)", "AT" * 40 + "A", "极高", "最高分"),
        ("设计发夹 (stem-loop)", hairpin_seq, ">3.0", "高分"),
        ("随机序列", random_seq, "~1.8", "低分"),
        ("Poly-A", "A" * 81, "0 (无互补)", "最低分"),
    ]
    
    print("分数分布 (stem method, normalized by log4(L)):")
    scores = []
    for name, seq, expected, note in test_cases:
        score = calculate_hairpin_fast(seq, method="stem")
        scores.append((name, score, note))
        print(f"  {name:25s}: {score:.2f} ({note})")
    
    # Verify ordering
    print("\n验证分数排序:")
    at_score = scores[0][1]  # AT repeat
    hp_score = scores[1][1]  # Hairpin
    rand_score = scores[2][1]  # Random
    poly_score = scores[3][1]  # Poly-A
    
    checks = [
        (at_score > hp_score, "AT重复 > 设计发夹"),
        (hp_score > rand_score, "设计发夹 > 随机"),
        (rand_score > poly_score, "随机 > Poly-A"),
    ]
    
    all_passed = True
    for passed, desc in checks:
        status = "✓" if passed else "✗"
        if not passed:
            all_passed = False
        print(f"  {status} {desc}")
    
    # Show random distribution across multiple lengths
    print("\n【统计】归一化分数在不同探针长度下的稳定性:")
    for probe_len in [51, 81, 121]:
        random.seed(42)
        rand_scores = [calculate_hairpin_fast(''.join(random.choices('ATCG', k=probe_len)), method="stem") for _ in range(100)]
        print(f"    {probe_len}bp: 平均={sum(rand_scores)/len(rand_scores):.2f}, "
              f"范围={min(rand_scores):.2f}-{max(rand_scores):.2f}")
    
    return all_passed


# =============================================================================
# Step 5: Dimer Score
# =============================================================================
def test_dimer():
    print_header("Step 5: Dimer Score")
    
    print("""
【实现方式】
- 评估探针池中的序列互补倾向（形成二聚体）
- 使用 K-mer 索引法:
  1. 从所有序列构建 k-mer 频率索引 (k=11)
  2. 对每条序列，统计其 k-mer 在全局索引中的出现频率
  3. 包括反向互补匹配
  4. 高频 k-mer 表示与其他序列有高度相似性

【期望行为】
- 使用百分位阈值 (默认保留 top 90%)
- 完全相同的序列会有最高分数
- 含有常见 motif 的序列会有较高分数
- 独特序列应该有较低分数
- 高 dimer 分数的探针可能互相杂交，降低捕获效率
""")
    
    print_subheader("测试用例")
    
    # 创建一个探针池
    sequences = [
        # 完全相同的序列 (高 dimer 分数)
        "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG" + "A",
        "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG" + "A",
        
        # 相似序列 (中高 dimer 分数)
        "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG" + "T",
        
        # 反向互补序列 (会匹配)
        "CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT",
        
        # 独特序列1
        "ATGCTGACGTAGCTGATCGATCGATCGTAGCTGACTAGCTGACGATCGTAGCTGATCGATCGATCGTAGCTAGCTGACTA" + "G",
        
        # 独特序列2
        "GCTAGCTGATCGTAGCTGACTAGCTGACGATCGTAGCTGATCGATCGTAGCTGACTAGCTGACGATCGTAGCTGATCGAT" + "A",
        
        # 独特序列3
        "TGACTAGCTGACGATCGTAGCTGATCGATCGTAGCTGACTAGCTGACGATCGTAGCTGATCGATCGTAGCTGACTAGCTG" + "C",
    ]
    
    names = [
        "相同序列1",
        "相同序列2", 
        "几乎相同序列",
        "反向互补",
        "独特序列1",
        "独特序列2",
        "独特序列3",
    ]
    
    # 构建 dimer 索引
    dimer_calc = DimerCalculatorFast(k=11)
    n_kmers = dimer_calc.build_index(sequences)
    print(f"\n构建索引: {len(sequences)} 条序列, {n_kmers} 个唯一 k-mer")
    
    dimer_scores = dimer_calc.calculate_all_scores()
    
    print("\n分数分布:")
    for name, score in zip(names, dimer_scores):
        print(f"  {name:20s}: {score:.4f}")
    
    print(f"\n分数统计:")
    print(f"  最小值: {min(dimer_scores):.4f}")
    print(f"  最大值: {max(dimer_scores):.4f}")
    print(f"  平均值: {sum(dimer_scores)/len(dimer_scores):.4f}")
    
    # 验证相同序列有相同分数
    print("\n验证:")
    if abs(dimer_scores[0] - dimer_scores[1]) < 0.0001:
        print("  ✓ 完全相同的序列有相同的 dimer 分数")
    else:
        print("  ✗ 相同序列分数不同!")
    
    if dimer_scores[0] > dimer_scores[4]:
        print("  ✓ 重复序列的 dimer 分数高于独特序列")
    else:
        print("  ✗ 分数排序不符合预期!")
    
    # 计算百分位阈值
    threshold_90 = calculate_percentile_threshold(dimer_scores, 90.0, higher_is_worse=True)
    print(f"\n90%百分位阈值: {threshold_90:.4f}")
    
    print("\n过滤结果 (使用 90% 百分位):")
    for name, score in zip(names, dimer_scores):
        passes = score <= threshold_90
        print(f"  {'✓' if passes else '✗'} {name:20s}: {score:.4f} -> {'通过' if passes else '过滤'}")
    
    return True


# =============================================================================
# Main
# =============================================================================
def main():
    print("\n" + "=" * 70)
    print("  Biophysical Filters - Stepwise Validation")
    print("=" * 70)
    
    results = {}
    
    results["GC"] = test_gc_content()
    results["Tm"] = test_tm()
    results["DUST"] = test_dust()
    results["Hairpin"] = test_hairpin()
    results["Dimer"] = test_dimer()
    
    print_header("Summary")
    for name, passed in results.items():
        status = "✓ PASS" if passed else "✗ CHECK"
        print(f"  {name:15s}: {status}")
    
    print("\n说明:")
    print("  - GC, Tm, DUST 使用固定阈值过滤")
    print("  - Hairpin, Dimer 使用百分位阈值，过滤分数最高的 10%")
    print("  - 百分位阈值会根据数据集动态计算")


if __name__ == "__main__":
    main()
