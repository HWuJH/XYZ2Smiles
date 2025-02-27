import streamlit as st
import pandas as pd
from rdkit import Chem

# 标题
st.title("XYZ 到 SMILES 转换器")
st.write("上传包含 XYZ 文件的文件夹，批量转换为 SMILES 并导出结果。")

# 上传文件夹
uploaded_files = st.file_uploader("上传 XYZ 文件", type="xyz", accept_multiple_files=True)

def extract_xyz_content(content):
    """
    从文件中提取标准的 XYZ 部分
    """
    lines = content.splitlines()
    # 找到原子坐标的开始行
    for i, line in enumerate(lines):
        if len(line.split()) == 4:  # 判断是否为原子坐标行
            xyz_lines = lines[i:]
            return "\n".join(xyz_lines)
    return None

if uploaded_files:
    st.write(f"已上传 {len(uploaded_files)} 个文件。")
    
    # 存储结果
    results = []
    
    # 遍历每个文件
    for uploaded_file in uploaded_files:
        try:
            # 读取文件内容
            content = uploaded_file.read().decode("utf-8")
            
            # 提取标准的 XYZ 部分
            xyz_content = extract_xyz_content(content)
            if not xyz_content:
                results.append({"文件名": uploaded_file.name, "SMILES": "文件格式错误"})
                continue
            
            # 使用 RDKit 将 XYZ 转换为分子对象
            mol = Chem.MolFromXYZBlock(xyz_content)
            if mol:
                # 将分子对象转换为 SMILES
                smiles = Chem.MolToSmiles(mol)
                results.append({"文件名": uploaded_file.name, "SMILES": smiles})
            else:
                results.append({"文件名": uploaded_file.name, "SMILES": "转换失败"})
        except Exception as e:
            results.append({"文件名": uploaded_file.name, "SMILES": f"错误: {str(e)}"})
    
    # 显示结果
    st.write("转换结果：")
    df = pd.DataFrame(results)
    st.dataframe(df)
    
    # 导出为 CSV
    if st.button("导出为 CSV"):
        csv = df.to_csv(index=False).encode("utf-8")
        st.download_button(
            label="下载 CSV",
            data=csv,
            file_name="smiles_results.csv",
            mime="text/csv",
        )
else:
    st.write("请上传 XYZ 文件。")
