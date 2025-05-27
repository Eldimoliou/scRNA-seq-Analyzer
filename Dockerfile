# Base image με Python
FROM python:3.10

# Ορίζουμε το working directory στο container
WORKDIR /app

# Αντιγράφουμε όλα τα αρχεία του project στο container
COPY . /app

# Εγκατάσταση των dependencies
RUN pip install --upgrade pip
RUN pip install -r requirements.txt

# Εκθέτουμε το port που χρησιμοποιεί το Streamlit
EXPOSE 8501

# Εκκίνηση της εφαρμογής
CMD ["streamlit", "run", "app.py", "--server.port=8501", "--server.address=0.0.0.0"]

COPY .streamlit /app/.streamlit
