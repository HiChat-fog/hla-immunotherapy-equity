import csv, os
os.makedirs("/home/mw/project/tables", exist_ok=True)

data = [
    ("A*02:01","European (NMDP)",27.1,"NMDP","欧洲裔美国人捐献者"),
    ("A*02:01","European (Romania)",26.1,"Caragea 2024","2024年高分辨率NGS数据"),
    ("A*02:01","European (median)",26.0,"AFND","跨欧洲多国中位频率"),
    ("A*02:01","East Asian (NMDP)",6.5,"NMDP","亚裔/太平洋岛民捐献者"),
    ("A*02:01","Chinese GC patients",5.0,"Zhou 2019","中国胃癌患者实测"),
    ("A*02:01","Japan",12.0,"AFND","日本人群"),
    ("A*02:01","China South Han",25.8,"AFND","浙江/广东汉族"),
    ("A*11:01","SE Asia (median)",21.0,"AFND","东南亚人群中位频率"),
    ("A*11:01","China South Han",27.7,"AFND","中国南方汉族"),
    ("A*11:01","Chinese GC patients",46.9,"Zhou 2019","中国胃癌患者实测"),
    ("A*11:01","China Yunnan Wa",58.4,"AFND","全球最高频率之一"),
    ("A*11:01","Japan",10.8,"AFND","日本频率"),
    ("A*11:01","South Korea",10.8,"AFND","韩国频率"),
    ("A*11:01","European (median)",5.6,"AFND","欧洲中位频率"),
    ("A*11:01","European (Romania)",8.6,"Caragea 2024","东南欧偏高"),
    ("A*11:01","Africa Sub-Saharan",0.1,"AFND","撒哈拉以南非洲几乎不存在"),
    ("A*24:02","Taiwan indigenous (Paiwan)",86.3,"AFND","南岛语系人群标志性等位基因"),
    ("A*24:02","Japan",32.7,"AFND","日本最常见HLA-A之一"),
    ("A*24:02","NE Asia (median)",23.0,"AFND","东北亚中位频率"),
    ("A*24:02","Chinese GC patients",25.0,"Zhou 2019","中国胃癌患者实测"),
    ("A*24:02","China North Han",15.2,"AFND","北汉频率"),
    ("A*24:02","European (median)",10.0,"AFND","欧洲中位频率"),
    ("A*24:02","European (Romania)",11.7,"Caragea 2024","巴尔干/东南欧"),
    ("A*24:02","NW Europe (Ireland)",6.8,"AFND","西北欧最低"),
    ("A*24:02","Czech Republic",1.0,"AFND","欧洲大陆最低点之一"),
]

out = "/home/mw/project/tables/hla_frequency_data.csv"
with open(out, "w", encoding="utf-8", newline="") as f:
    w = csv.writer(f)
    w.writerow(["allele","population","frequency_pct","source","notes"])
    w.writerows(data)
print(f"已生成: {out} ({len(data)} 条)")
