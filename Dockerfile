FROM python:3.11.9-slim-bookworm

LABEL maintainer="Juniper Lake <jlake@pacb.com>"

ARG PACKAGE="humanatee"

COPY --from=ghcr.io/astral-sh/uv:latest /uv /bin/uv

WORKDIR /app

RUN --mount=type=bind,source=uv.lock,target=uv.lock \
    --mount=type=bind,source=pyproject.toml,target=pyproject.toml \
    uv sync --frozen --no-install-project --no-dev --package=$PACKAGE

COPY src/$PACKAGE /app/src/$PACKAGE

RUN --mount=type=bind,source=uv.lock,target=uv.lock \
    --mount=type=bind,source=pyproject.toml,target=pyproject.toml \
    uv sync --frozen --no-dev --package=$PACKAGE

ENV PATH="/app/.venv/bin:$PATH"

RUN humanatee version
