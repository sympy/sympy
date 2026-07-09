import asyncio
import json
import subprocess
import signal
import os

from contextlib import suppress

import httpx
from fastapi import FastAPI, Request
from fastapi.responses import JSONResponse, StreamingResponse
from fastapi.middleware.cors import CORSMiddleware

app = FastAPI()

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# -------------------------------------------------------
# MODEL PROFILES
# -------------------------------------------------------

MODEL_PROFILES = {
    "gemma4-31b-qat": {
        "path": "/home/grdi/.lmstudio/models/unsloth/gemma-4-31B-it-qat-GGUF/gemma-4-31B-it-qat-UD-Q4_K_M.gguf",
        "backend": "llama",  # "llama" or "ik_llama"
        "ctx": "32768",
        "layers": "99",
        "cache_k": None,
        "cache_v": None,
        "type": "chat",
        "server_args": {
            "flash-attn": "on",
            "temp": 1.1,
            "no-kv-offload": False,
        },
    },

    "gemma4-31b-qat-kv_q8-q6": {
        "path": "/home/grdi/.lmstudio/models/unsloth/gemma-4-31B-it-qat-GGUF/gemma-4-31B-it-qat-UD-Q4_K_M.gguf",
        "backend": "llama",
        "ctx": "32768",
        "layers": "99",
        "cache_k": "q6_0",
        "cache_v": "q8_0",
        "type": "chat",
        "server_args": {
            "flash-attn": "on",
            "temp": 1.1,
            "no-kv-offload": True,
        },
    },

    "qwen3.6-27b-q4_k_xl-hightemp": {
        "path": "/home/grdi/.lmstudio/models/unsloth/Qwen3.6-27B-GGUF/Qwen3.6-27B-UD-Q4_K_XL.gguf",
        "backend": "ik_llama",
        "ctx": "99000",
        "layers": "99",
        "cache_k": None,
        "cache_v": None,
        "type": "chat",
        "server_args": {
            "flash-attn": "on",
            "temp": 1,
            "no-kv-offload": False,
            "metrics": True,
            "jinja": True,
        },
    },

    "qwen3.6-27b-q4_k_xl-lowtemp": {
        "path": "/home/grdi/.lmstudio/models/unsloth/Qwen3.6-27B-GGUF/Qwen3.6-27B-UD-Q4_K_XL.gguf",
        "backend": "ik_llama",
        "ctx": "99000",
        "layers": "99",
        "cache_k": None,
        "cache_v": None,
        "type": "chat",
        "server_args": {
            "flash-attn": "on",
            "temp": 0.6,
            "no-kv-offload": False,
            "metrics": True,
            "jinja": True,
        },
    },

    "qwen3.6-27b-q6_k-hightemp": {
        "path": "/home/grdi/.lmstudio/models/unsloth/Qwen3.6-27B-GGUF/Qwen3.6-27B-Q6_K.gguf",
        "backend": "ik_llama",
        "ctx": "150000",
        "layers": "99",
        "cache_k": None,
        "cache_v": None,
        "type": "chat",
        "server_args": {
            "flash-attn": "on",
            "temp": 1,
            "no-kv-offload": True,
            "metrics": True,
            "jinja": True,
        },
    },

    "qwen3.6-27b-q6_k-lowtemp": {
        "path": "/home/grdi/.lmstudio/models/unsloth/Qwen3.6-27B-GGUF/Qwen3.6-27B-Q6_K.gguf",
        "backend": "ik_llama",
        "ctx": "150000",
        "layers": "99",
        "cache_k": None,
        "cache_v": None,
        "type": "chat",
        "server_args": {
            "flash-attn": "on",
            "temp": 0.6,
            "no-kv-offload": True,
            "metrics": True,
            "jinja" : True,
        },
    },

     "qwen3.5-122b-lowtemp": {
        "path": "/home/grdi/.lmstudio/models/unsloth/Qwen3.5-122B-A10B-GGUF/Qwen3.5-122B-A10B-Q8_0-00001-of-00004.gguf",
        "backend": "ik_llama",
        "ctx": "32768",
        "layers": "1",
        "cache_k": None,
        "cache_v": None,
        "type": "chat",
        "no_cuda_pinned_memory": False,

        "server_args": {
            "override-tensor": "exps=CPU",
            "flash-attn": "on",
            "temp": 0.6,
            "batch-size": "256",
            "ubatch-size": "128",
            "metrics": True,
            "jinja": True,
        },
    },

    "qwen3.5-122b-hightemp": {
        "path": "/home/grdi/.lmstudio/models/unsloth/Qwen3.5-122B-A10B-GGUF/Qwen3.5-122B-A10B-Q8_0-00001-of-00004.gguf",
        "backend": "ik_llama",
        "ctx": "450000",
        "layers": "99",
        "cache_k": None,
        "cache_v": None,
        "type": "chat",
        "server_args": {
            "flash-attn": "on",
            "temp": 1,
            "batch-size": "8192",
            "ubatch-size": "2048",
            "cpu-moe": True,
            "override-tensor": "exps=CPU",
            "metrics": True,
            "jinja": True,
        },
    },

    "qwen3-embedding-8b": {
        "path": "/home/grdi/.lmstudio/models/Qwen/Qwen3-Embedding-8B-GGUF/Qwen3-Embedding-8B-f16.gguf",
        "backend": "llama",
        "ctx": "32768",
        "layers": "99",
        "cache_k": None,
        "cache_v": None,
        "type": "embedding",
    },
}

# Absolute directory paths relative to your layout
BACKEND_BINARIES = {
    "llama": "/home/grdi/ai/tools/llama.cpp/build/bin/llama-server",
    "ik_llama": "/home/grdi/ai/tools/ik_llama.cpp/build/bin/llama-server"
}

INTERNAL_PORT = 8080
PROXY_PORT = 1234

http_client = httpx.AsyncClient(
    base_url=f"http://127.0.0.1:{INTERNAL_PORT}",
    timeout=None,
)

server_process = None
current_config_hash = None
reload_lock = asyncio.Lock()


# -------------------------------------------------------
# SERVER CONTROL
# -------------------------------------------------------

async def stop_server():
    global server_process

    if not server_process:
        return

    print("Stopping llama-server...")

    with suppress(Exception):
        server_process.send_signal(signal.SIGINT)
        try:
            server_process.wait(timeout=10)
        except Exception:
            server_process.terminate()
            try:
                server_process.wait(timeout=5)
            except Exception:
                server_process.kill()

    server_process = None


async def start_llama_server(profile_name, ctx_override=None, layers_override=None,
                             cache_k_override=None, cache_v_override=None):

    global server_process, current_config_hash

    async with reload_lock:

        config = MODEL_PROFILES.get(profile_name)
        if not config:
            raise ValueError(f"Unknown profile: {profile_name}")

        backend_type = config.get("backend", "llama")
        binary_path = BACKEND_BINARIES.get(backend_type)
        if not binary_path:
            raise ValueError(f"Unknown backend type '{backend_type}' specified for profile: {profile_name}")

        is_embedding = config.get("type") == "embedding"

        # Override values only if explicitly provided (handle 0 correctly)
        final_ctx = str(ctx_override) if ctx_override is not None else config["ctx"]
        final_kv = str(layers_override) if layers_override is not None else config["layers"]
        final_k = cache_k_override if cache_k_override is not None else config.get("cache_k")
        final_v = cache_v_override if cache_v_override is not None else config.get("cache_v")

        # Effective GPU layers (embedding always uses 0)
        gpu_layers = "99" if is_embedding else final_kv

        # Hash based on the *effective* configuration and chosen backend repository
        config_hash = hash((profile_name, backend_type, final_ctx, gpu_layers, final_k, final_v))

        if server_process and config_hash == current_config_hash:
            return

        await stop_server()

        print(f"Starting {profile_name} using backend: {backend_type} (embedding={is_embedding})")

        cmd = [
            binary_path,
            "--model", config["path"],
            "--host", "0.0.0.0",
            "--port", str(INTERNAL_PORT),
            "--ctx-size", final_ctx,
            "--n-gpu-layers", gpu_layers,
            "--parallel", "1",
            "--threads", "6",
        ]

        if is_embedding:
            cmd.append("--embedding")
        if final_k:
            cmd.extend(["--cache-type-k", final_k])
        if final_v:
            cmd.extend(["--cache-type-v", final_v])

        # Generic server arguments
        if not is_embedding:
            for key, value in config.get("server_args", {}).items():
                flag = f"--{key}"

                if isinstance(value, bool):
                    if value:
                        cmd.append(flag)
                elif value is not None:
                    cmd.extend([flag, str(value)])

        print("CMD:", " ".join(cmd))
        
        env = dict(os.environ)
        
        if config.get("no_cuda_pinned_memory", False):
            env["GGML_CUDA_NO_PINNED"] = "1"
        
        server_process = subprocess.Popen(
            cmd,
            env = env,
            )
            
        current_config_hash = config_hash

        print("Waiting for llama-server...")

        # Wait for server to become ready, but with a timeout
        max_wait = 120  # seconds
        waited = 0
        while True:
            poll = server_process.poll()
            if poll is not None:
                raise RuntimeError(f"llama-server crashed immediately (exit {poll})")

            try:
                r = await http_client.get("/v1/models")
                if r.status_code == 200:
                    break
            except Exception:
                pass

            await asyncio.sleep(0.5)
            waited += 0.5
            if waited >= max_wait:
                server_process.kill()
                server_process = None
                raise TimeoutError("llama-server did not start within 120 seconds")

        print("Ready.")


# -------------------------------------------------------
# PROXY ENDPOINTS
# -------------------------------------------------------

@app.get("/v1/models")
async def list_models():
    return {
        "object": "list",
        "data": [
            {
                "id": name,
                "object": "model",
                "created": 1700000000,
                "owned_by": "custom-manager"
            }
            for name in MODEL_PROFILES.keys()
        ]
    }


@app.api_route("/{path:path}", methods=["GET", "POST", "PUT", "DELETE"])
async def proxy(request: Request, path: str):

    try:
        # Read the body once – prevents double‑consumption bug
        body = await request.body()

        # If it's a POST with a JSON body, parse it to check for a model change
        if request.method == "POST" and body:
            try:
                body_json = json.loads(body)
            except Exception:
                body_json = None

            if body_json and "model" in body_json:
                model = body_json["model"]

                if model not in MODEL_PROFILES:
                    return JSONResponse({"error": f"Unknown model {model}"}, status_code=400)

                await start_llama_server(
                    model,
                    body_json.get("ctx_override"),
                    body_json.get("layers_override"),
                    body_json.get("cache_k_override"),
                    body_json.get("cache_v_override"),
                )

        # Forward the request to the internal llama‑server
        headers = dict(request.headers)
        headers.pop("host", None)
        headers.pop("content-length", None)
        headers.pop("connection", None)

        url = httpx.URL(path="/" + path, query=request.url.query.encode("utf-8"))

        req = http_client.build_request(
            method=request.method,
            url=url,
            headers=headers,
            content=body,
        )

        resp = await http_client.send(req, stream=True)

        response_headers = dict(resp.headers)
        response_headers.pop("content-length", None)
        response_headers.pop("transfer-encoding", None)

        return StreamingResponse(
            resp.aiter_raw(),
            status_code=resp.status_code,
            headers=response_headers,
            media_type=resp.headers.get("content-type"),
        )

    except Exception as e:
        return JSONResponse({"error": str(e)}, status_code=502)


# -------------------------------------------------------
# SHUTDOWN
# -------------------------------------------------------

@app.on_event("shutdown")
async def shutdown():
    await http_client.aclose()
    await stop_server()


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=PROXY_PORT)

