.PHONY: install
install:
	uv run fpm install --profile release --prefix .
