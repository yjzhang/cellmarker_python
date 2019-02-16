nohup gunicorn --workers 4 --threads 1 --bind 0.0.0.0:8889 --log-level debug -t 5001 cellmarker_app:app > nohup_cellmarker.out &

