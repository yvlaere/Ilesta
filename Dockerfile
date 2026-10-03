# Stage 1: Build Ilesta and fetch minimap2 via Pixi
FROM ghcr.io/prefix-dev/pixi:0.81.0 AS build

WORKDIR /app

# Copy dependency specifications and install environments
COPY pixi.toml pixi.lock ./
RUN pixi install --locked -e default && \
    pixi install --locked -e dev

# Copy source code and build Ilesta release binary
COPY Cargo.toml Cargo.lock ./
COPY src/ ./src/
RUN pixi run -e dev cargo build --release

# Stage 2: Minimal production runtime image (no entrypoints, ~80 MB total)
FROM debian:bookworm-slim AS production

# Copy compiled Ilesta binary and pixi-managed minimap2 binary directly into /usr/local/bin
COPY --from=build /app/target/release/Ilesta /usr/local/bin/Ilesta
COPY --from=build /app/.pixi/envs/default/bin/minimap2 /usr/local/bin/minimap2

# Copy minimap2's dynamic shared library dependency (libz) and register with linker
COPY --from=build /app/.pixi/envs/default/lib/libz.so.1* /usr/local/lib/
RUN ldconfig

WORKDIR /data

CMD ["Ilesta", "--help"]
