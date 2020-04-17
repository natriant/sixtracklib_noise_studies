#!/bin/bash

# run the job files
for f in job*.py; do
  echo running "$f"
  python "$f" `echo $@` || break  # execute successfully or break
done