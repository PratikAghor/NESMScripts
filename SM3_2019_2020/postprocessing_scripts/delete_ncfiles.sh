#!/bin/bash
for dir in aghor_*; do
  rm -rf "$dir"/*.nc "$dir"/Report-*
done
