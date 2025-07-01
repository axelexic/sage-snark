#!/bin/sh

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

VENV_DIR_NAME="$SCRIPT_DIR/venv"
VENV_PROMPT_NAME=$(basename $SCRIPT_DIR)

SAGE="$HOME/Applications/SageMath-10-6.app/Contents/Frameworks/Sage.framework/Versions/10.6/local/var/lib/sage/venv-python3.12.5/bin/sage"

if [ -d "$VENV_DIR_NAME" ]; then
  echo "Virtual environment \"$VENV_DIR_NAME\" already exists"
else
	set -x
	"$SAGE" --python -m venv --system-site-packages --prompt "$VENV_PROMPT_NAME" "$VENV_DIR_NAME"
fi
