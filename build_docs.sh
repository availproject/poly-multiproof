#!/bin/bash
RUSTDOCFLAGS="--html-in-header doc/katex-header.html" cargo doc --document-private-items
