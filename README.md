## Troubleshoot

If `cargo build` is not working, and the compiler indicates that there is a missing library (`/usr/bin/ld: cannot find -lopenblas: No such file or directory`), try install the library with

```bash
sudo apt install libopenblas-dev
```