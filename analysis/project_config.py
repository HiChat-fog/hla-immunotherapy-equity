"""
============================================================
  项目共享配置
  国家名映射、区域划分、驱动基因列表等
  所有步骤脚本共用此文件
============================================================
"""

# ====== 中文国家名 → 英文（匹配 World Bank 数据） ======
CN_TO_EN = {
    # ---- 东亚与太平洋 ----
    "中国": "China",
    "日本": "Japan",
    "韩国": "Korea, Rep.",
    "大韩民国": "Korea, Rep.",
    "越南": "Viet Nam",
    "泰国": "Thailand",
    "菲律宾": "Philippines",
    "印度尼西亚": "Indonesia",
    "马来西亚": "Malaysia",
    "蒙古": "Mongolia",
    "缅甸": "Myanmar",
    "柬埔寨": "Cambodia",
    "老挝": "Lao PDR",
    "老挝人民民主共和国": "Lao PDR",
    "新加坡": "Singapore",
    "斐济": "Fiji",
    # ---- 欧洲与中亚 ----
    "德国": "Germany",
    "法国": "France",
    "英国": "United Kingdom",
    "意大利": "Italy",
    "西班牙": "Spain",
    "俄罗斯": "Russian Federation",
    "俄罗斯联邦": "Russian Federation",
    "波兰": "Poland",
    "荷兰": "Netherlands",
    "比利时": "Belgium",
    "瑞典": "Sweden",
    "瑞士": "Switzerland",
    "奥地利": "Austria",
    "希腊": "Greece",
    "葡萄牙": "Portugal",
    "丹麦": "Denmark",
    "芬兰": "Finland",
    "挪威": "Norway",
    "爱尔兰": "Ireland",
    "捷克共和国": "Czechia",
    "匈牙利": "Hungary",
    "罗马尼亚": "Romania",
    "乌克兰": "Ukraine",
    "保加利亚": "Bulgaria",
    "塞尔维亚": "Serbia",
    "克罗地亚": "Croatia",
    "斯洛伐克": "Slovak Republic",
    "立陶宛": "Lithuania",
    "白俄罗斯": "Belarus",
    "哈萨克斯坦": "Kazakhstan",
    "阿塞拜疆": "Azerbaijan",
    "亚美尼亚": "Armenia",
    "摩尔多瓦共和国": "Moldova",
    # ---- 北美 ----
    "美利坚合众国": "United States",
    "美国": "United States",
    "加拿大": "Canada",
    # ---- 拉丁美洲 ----
    "巴西": "Brazil",
    "墨西哥": "Mexico",
    "阿根廷": "Argentina",
    "哥伦比亚": "Colombia",
    "智利": "Chile",
    "秘鲁": "Peru",
    "古巴": "Cuba",
    "委内瑞拉玻利瓦尔共和国": "Venezuela, RB",
    # ---- 南亚 ----
    "印度": "India",
    "巴基斯坦": "Pakistan",
    "孟加拉国": "Bangladesh",
    "尼泊尔": "Nepal",
    "斯里兰卡": "Sri Lanka",
    # ---- 撒哈拉以南非洲 ----
    "尼日利亚": "Nigeria",
    "南非": "South Africa",
    "肯尼亚": "Kenya",
    "埃塞俄比亚": "Ethiopia",
    "坦桑尼亚": "Tanzania",
    "坦桑尼亚联合共和国": "Tanzania",
    "刚果民主共和国": "Congo, Dem. Rep.",
    "乌干达": "Uganda",
    "加纳": "Ghana",
    "莫桑比克": "Mozambique",
    "安哥拉": "Angola",
    "喀麦隆": "Cameroon",
    "科特迪瓦": "Cote d'Ivoire",
    # ---- 中东与北非 ----
    "沙特阿拉伯": "Saudi Arabia",
    "伊朗伊斯兰共和国": "Iran, Islamic Rep.",
    "埃及阿拉伯共和国": "Egypt, Arab Rep.",
    "伊拉克": "Iraq",
    "摩洛哥": "Morocco",
    "阿尔及利亚": "Algeria",
    "阿拉伯联合酋长国": "United Arab Emirates",
    "卡塔尔": "Qatar",
    "科威特": "Kuwait",
    "土耳其": "Turkiye",
    # ---- 大洋洲 ----
    "澳大利亚": "Australia",
    "新西兰": "New Zealand",
    # ---- 其他 ----
    "苏丹": "Sudan",
    "朝鲜民主主义人民共和国": "Korea, Dem. People's Rep.",
}

# ====== 区域 → 中文国家名列表 ======
REGIONS = {
    "东亚与太平洋": [
        "中国", "日本", "韩国", "大韩民国", "越南", "泰国", "菲律宾",
        "印度尼西亚", "马来西亚", "蒙古", "缅甸", "柬埔寨",
        "老挝", "老挝人民民主共和国", "斐济", "新加坡"
    ],
    "欧洲与中亚": [
        "德国", "法国", "英国", "意大利", "西班牙", "俄罗斯", "俄罗斯联邦",
        "波兰", "荷兰", "比利时", "瑞典", "瑞士", "奥地利", "希腊", "葡萄牙",
        "丹麦", "芬兰", "挪威", "爱尔兰", "捷克共和国", "匈牙利", "罗马尼亚",
        "乌克兰", "保加利亚", "塞尔维亚", "克罗地亚", "斯洛伐克", "立陶宛"
    ],
    "北美": ["美利坚合众国", "美国", "加拿大"],
    "拉丁美洲与加勒比": [
        "巴西", "墨西哥", "阿根廷", "哥伦比亚", "智利", "秘鲁",
        "古巴", "委内瑞拉玻利瓦尔共和国"
    ],
    "南亚": ["印度", "巴基斯坦", "孟加拉国", "尼泊尔", "斯里兰卡"],
    "撒哈拉以南非洲": [
        "尼日利亚", "南非", "肯尼亚", "埃塞俄比亚", "坦桑尼亚", "坦桑尼亚联合共和国",
        "刚果民主共和国", "乌干达", "加纳", "莫桑比克", "安哥拉", "喀麦隆", "科特迪瓦"
    ],
    "中东与北非": [
        "沙特阿拉伯", "伊朗伊斯兰共和国", "埃及阿拉伯共和国", "伊拉克",
        "摩洛哥", "阿尔及利亚", "阿拉伯联合酋长国", "卡塔尔", "科威特"
    ],
}

# ====== 胃癌驱动基因 ======
DRIVER_GENES = {
    "TP53", "ARID1A", "KRAS", "PIK3CA", "CDH1", "ERBB2", "RHOA",
    "SMAD4", "CTNNB1", "APC", "PTEN", "FBXW7", "ERBB3", "MET"
}

# ====== 输出目录 ======
import os
FIG_DIR = os.path.join(os.path.dirname(__file__), "figures")
