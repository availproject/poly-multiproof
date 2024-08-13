#!/bin/bash
RUSTDOCFLAGS="--html-in-header $(pwd)/doc/katex-header.html" cargo doc --document-private-items
