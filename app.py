import streamlit as st
import pandas as pd
from openbabel import pybel

# 标题
st.title("XYZ 到 SMILES 转换器")
st.write("上传包含 XYZ 文件的文件夹，批量转换为 SMILES 并导出结果。")

# 上传文件夹
uploaded_files = st.file_uploader("上传 XYZ 文件", type="xyz", accept_multiple_files=True)

if uploaded_files:
    st.write(f"已上传 {len(uploaded_files)} 个文件。")
    
    # 存储结果
    results = []
    
    # 遍历每个文件
    for uploaded_file in uploaded_files:
        try:
            # 读取文件内容
            content = uploaded_file.read().decode("utf-8")
            
            # 使用 Open Babel 将 XYZ 转换为分子对象
            mol = pybel.readstring("xyz", content)
            if mol:
                # 将分子对象转换为 SMILES
                smiles = mol.write("smi")
                results.append({"文件名": uploaded_file.name, "SMILES": smiles.strip()})
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
