# PICOTA release checklist

The steps for cutting a release. Version history and what changed in each
release are in [CHANGELOG.md](CHANGELOG.md).

## 1. Before tagging

- [ ] The full suite passes, with the BLAST databases in place, so the
      end-to-end tests run instead of skipping:

      ```bash
      conda run -n evobiomig python -m pytest -q
      ```

      The `xfail` entries are the open defects in `tests/test_known_defects.py`
      and are expected. An `XPASS` failure means a defect was fixed and its
      marker has to be removed.

- [ ] The version is the same in all three places, which `tests/test_cli.py`
      checks for the first two:
  - `setup.py`
  - `picota/picota.py` (`__version__` and the docstring)
  - `picota/__init__.py`

- [ ] `CHANGELOG.md` has a dated section for the release, and every change that
      alters what an existing run writes is listed under **Output changes**.

- [ ] Every path in `picota/config.yaml` exists, and the README's quick start
      runs as written from the repository root:

      ```bash
      python3 picota/picota.py analysis --gfa picota/test_data/testNitro.gfa --output /tmp/picota_check
      ```

- [ ] If a reference database or the scoring changed, regenerate the benchmark
      figures quoted in the README with `scripts/build_report.sh`, or mark them
      as predating the change.

- [ ] The package installs and its entry point imports:

      ```bash
      pip install --no-deps --target /tmp/picota_inst . && PYTHONPATH=/tmp/picota_inst /tmp/picota_inst/bin/picota --version
      ```

## 2. Tag and push

The remote uses SSH (`git@github.com:recepcanaltinbag/picota.git`).

```bash
git push origin main
git tag -a vX.Y.Z -m "PICOTA vX.Y.Z"
git push origin vX.Y.Z
```

## 3. GitHub release

1. <https://github.com/recepcanaltinbag/picota/releases> → *Draft a new release*.
2. Select the tag and title the release `PICOTA vX.Y.Z`.
3. Paste the release's CHANGELOG section as the description, with its
   **Output changes** first.
4. The reference databases are not in the repository. If they are distributed
   with the release, check each source's terms (CARD, ISfinder, TnCentral,
   KEGG, UniProt) before attaching them.
5. If the Zenodo record in the README is linked to GitHub releases, check that
   the new version appears there. Otherwise upload it to Zenodo by hand.
