import time
import os
import re
import glob

try:
    import gmpy2
    HAS_GMPY2 = True
except ImportError:
    HAS_GMPY2 = False
    print("[!] 警告：请安装 gmpy2 以获得极速体验 (pip install gmpy2)")

Q = 47
# 扫描主星左侧多远的距离？(步长为2的偶数，建议 5000)
SEARCH_RADIUS = 5000

def q47(n):
    return n**Q - (n-1)**Q

def main():
    print("==================================================")
    print("📡 泰坦深空伴星雷达 v3.0 (终极全星表阵列)")
    print(f"🎯 扫描半径: P - 2 到 P - {SEARCH_RADIUS}")
    print("==================================================")

    n_values = set()
    
    # 1. 自动读取当前目录下所有的日志文件 (捕捉深空新发现)
    all_files = glob.glob("*.log") + glob.glob("*.txt")
    
    for filename in all_files:
        try:
            with open(filename, 'r', encoding='utf-8') as f:
                content = f.read()
                matches = re.findall(r'(?:Sequence|QUADRUPLET):\s*(\d+)', content)
                for m in matches:
                    n_values.add(int(m))
        except Exception as e:
            pass

    # 2. 注入核心数据：0-40 亿硬编码星表
    hardcoded_n = [
        # --- 附件表格中的 14 颗元老级 4连星 (0 - 20亿) ---
        117309848, 136584738, 218787064, 411784485, 423600750, 
        523331634, 640399031, 987980498, 1163461515, 1370439187, 
        1643105964, 1691581855, 1975860550, 1996430175,
        
        # --- 聊天记录中的 11 颗中期 4连星 (20亿 - 40亿) ---
        2156109985, 2367719045, 2559344807, 2646631730, 2682956949, 
        2859276863, 2862155914, 2922108368, 3808591354, 3910149357, 
        3984049296
    ]
    
    for n in hardcoded_n:
        n_values.add(n)

    sorted_n = sorted(list(n_values))
    total_main_stars = len(sorted_n) * 4
    
    print(f"[!] 成功点亮 {len(sorted_n)} 组巨型雷达基站")
    print(f"[!] 即将对 {total_main_stars} 颗极其罕见的主星周边发起深空探测...")
    print("--------------------------------------------------")

    if not HAS_GMPY2:
        print("[!] 警告：未检测到 gmpy2，高维素性测试将会非常缓慢！")

    total_satellites = 0
    twin_primes = 0

    # 3. 开始雷达扫描
    start_time = time.time()
    
    for base_n in sorted_n:
        for offset in range(4): # 遍历 4连星 的每一颗主星
            n = base_n + offset
            P = q47(n)
            
            # 向左侧撒网
            for k in range(2, SEARCH_RADIUS + 1, 2):
                # 【泰坦护盾过滤器】：直接跳过必定被 3 整除的死区
                if k % 3 == 1:
                    continue
                
                candidate = P - k
                
                # 极速素性测试
                if HAS_GMPY2 and gmpy2.is_prime(candidate, 25):
                    total_satellites += 1
                    if k == 2:
                        twin_primes += 1
                        print(f"🚨🚨 [世纪发现！] 捕获纯血孪生素数！主星 n = {n}, 伴星 = P - 2")
                    else:
                        print(f"🛰️  发现伴星！主星 n = {n}, 伴星间距 = P - {k}")
                        
    end_time = time.time()
    
    print("==================================================")
    print(f"🏁 雷达扫描完毕！耗时: {end_time - start_time:.2f} 秒")
    print(f"共勘测了 {total_main_stars} 颗主星附近的引力空域。")
    print(f"🎯 总计捕获卫星：{total_satellites} 颗！")
    print(f"👑 纯血孪生素数 (P, P-2)：{twin_primes} 对！")

if __name__ == "__main__":
    main()