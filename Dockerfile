# Serves already-converted data (data/structures/*.cif + data/index/*.parquet).
# Data conversion (scripts/build_data.py) must be run beforehand on the host,
# where bin/foldcomp -- a macOS build -- can actually run; the container never
# calls foldcomp itself.

FROM python:3.11-slim

WORKDIR /app

COPY pyproject.toml requirements.txt ./
COPY src/ src/
RUN pip install --no-cache-dir -r requirements.txt \
    && pip install --no-cache-dir -e .

COPY app.py ./
COPY pages/ pages/
COPY .streamlit/ .streamlit/

# Mount the converted data/ directory here at run time, e.g.:
#   docker run -p 8501:8501 -v /path/to/data:/data eukaryoma-ppi-webserver
ENV EUKARYOMA_DATA_DIR=/data

EXPOSE 8501

ENTRYPOINT ["streamlit", "run", "app.py", "--server.port=8501", "--server.address=0.0.0.0"]
