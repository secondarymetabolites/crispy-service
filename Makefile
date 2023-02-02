unit:
	pytest -v

coverage:
	pytest --cov=crispy_service --cov-report=html --cov-report=term-missing


lint:
	flake8 crispy_service --count --select=E9,F63,F7,F82 --show-source --statistics
	flake8 crispy_service --count --exit-zero --max-complexity=10 --max-line-length=127 --statistics
	mypy crispy_service --ignore-missing-imports
