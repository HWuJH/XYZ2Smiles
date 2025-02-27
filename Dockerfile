# 使用 Python 3.9 作为基础镜像
FROM python:3.9-slim

# 安装系统依赖（包括 Open Babel）
RUN apt-get update && apt-get install -y \
    openbabel \
    && rm -rf /var/lib/apt/lists/*

# 安装 Python 依赖
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# 复制应用代码
COPY . .

# 运行 Streamlit 应用
CMD ["streamlit", "run", "app.py", "--server.port=8501", "--server.address=0.0.0.0"]
