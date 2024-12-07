import re  
import csv  

# Tên tệp đầu vào và đầu ra  
input_file = "asia_data/data/nextstrain.nh"  # File chứa dữ liệu Newick  
output_file = "asia_data/asia_epi_isl_ids_sorted.csv"  # File CSV lưu kết quả  

# Đọc dữ liệu từ tệp newick.nh  
with open(input_file, "r") as file:  
    data = file.read()  

# Biểu thức chính quy để tìm tất cả các ID EPI_ISL_  
pattern = r"EPI_ISL_\d+"  
epi_isl_ids = re.findall(pattern, data)  

# Loại bỏ các ID trùng lặp (nếu có) và sắp xếp theo thứ tự tăng dần  
unique_sorted_ids = sorted(set(epi_isl_ids), key=lambda x: int(x.split('_')[-1]))  

# Ghi kết quả vào file CSV  
with open(output_file, mode="w", newline="") as csvfile:  
    writer = csv.writer(csvfile)  
    writer.writerow(["EPI_ISL_ID"])  # Ghi tiêu đề cột  
    for epi_id in unique_sorted_ids:  
        writer.writerow([epi_id])  

print(f"Kết quả đã được lưu vào tệp '{output_file}'.") 