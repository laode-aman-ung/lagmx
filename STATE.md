# Status

Diperbarui: 2026-09-02 — PC Ubuntu 24 (`arga`)

## Sedang dikerjakan

Penggantian nama `pyAutoGMX` → **LAGMX** selesai dan terdorong. Repo bersih.

## Yang berubah hari ini

- Nama diganti di seluruh proyek: 82 kemunculan pada 16 berkas, kedua ejaan
  (`pyAutoGMX` dan `pyautogmx`). Berkas ikut di-rename: `LAGMX.py`,
  `run_matrix/run_lagmx.sh`.
- Judul naskah JOSS dan `CITATION.cff` (title, repository-code, url) mengikuti.
- Repo GitHub di-rename ke `laode-aman-ung/lagmx`; GitHub mengalihkan URL lama
  sehingga tautan yang sudah tersebar tetap hidup.
- Folder lokal pindah ke `~/riset/LAGMX`, mengikuti pola path identik lintas
  mesin.
- Ditambahkan berkas wajib sinkronisasi: `CLAUDE.md`, `STATE.md`,
  `.gitattributes`.

## JOSS — diajukan ulang 2 September 2026

**Status: submitted, review belum dimulai.**
https://joss.theoj.org/papers/b5d7f43007f1198a2e40430e524621a1

Diajukan sebagai *Resubmission* dari #11216, versi **v1.1.0**, subjek
Computational chemistry. Catatan untuk editor (2.973 karakter) memuat keenam
tautan repo pendahulu, argumen bahwa tanggal pembuatan repo adalah metadata
sisi server, tabel perbandingan urutan perintah `gmx` antara skrip 2023 dan
`LAGMX.py`, serta pengakuan terbuka atas kelemahannya: keenam repo diunggah
sekali, ada jeda 2023–2026, dan LAGMX adalah tulis ulang.

**Menunggu triase.** Email konfirmasi JOSS menyatakan naskah menunggu triase
managing editor — *"generally takes up to a week, but during busy times can
take a number of weeks"* — dan sampai lolos desk check, naskah **tidak muncul**
di `openjournals/joss-reviews`. Jadi tidak ada yang bisa dipantau; pemberitahuan
datang lewat email `[PRE REVIEW]` dari GitHub.

**Begitu issue muncul**, tempel komentar dari
`~/riset/_arsip/lagmx-joss-review-comment.md`. Itu bukan sekadar merapikan
format: isi *Notes to editor* tidak ikut tampil di badan issue publik — pada
#11216 issue-nya hanya memuat metadata — sehingga seluruh argumen silsilah
hanya terlihat di panel admin editor sampai komentar itu ditempel.

Bila ditolak lagi, tanggal paling cepat berikutnya tetap **24 Februari 2027**
(enam bulan sejak commit pertama 26 Agustus 2026).

## Riwayat penolakan pertama — 27 Agustus 2026

Pengajuan sebagai pyAutoGMX (issue #11216) **ditolak oleh AEiC Rachel
Kurchin**, bukan karena galat teknis:

> "I'll need to reject this submission under our new scope guidelines. Your
> code needs to have been under open development for at least six months."

Seluruh pemeriksaan otomatis justru lolos: 7 DOI OK tanpa satu pun MISSING
atau INVALID, lisensi MIT diakui OSI, Statement of need ditemukan, proof PDF
terbentuk. Yang menjatuhkan adalah umur dan sinyal aktivitas:

```
commit pertama 26 Agustus 2026 · diajukan 27 Agustus 2026  → berumur 1 hari
896 insertions dalam satu jendela 48 jam                    ditandai merah
stars 0 · forks 0 · contributors 1 · releases 0 · issues 0 · PR 0
```

**Tanggal paling cepat: 24 Februari 2027** (enam bulan sejak commit pertama).

Penggantian nama **tidak mereset jam itu** — riwayat git utuh dan tanggal
pembuatan repo di GitHub tetap 26 Agustus 2026. Membuat repo baru alih-alih
me-rename akan memulai hitungan dari nol; itu sebabnya rename dipilih.

### Yang harus dibangun sampai Februari

Enam bulan adalah syarat minimum, bukan jaminan. Yang dinilai editor adalah
apakah ini perangkat lunak yang hidup, dan semua angka di atas nol.

1. **Commit tersebar**, bukan menumpuk. Pola satu ledakan lalu senyap persis
   yang ditandai merah oleh editorialbot.
2. **Rilis bertahap.** Mulai dari versi yang menandai penggantian nama, lalu
   seterusnya mengikuti perbaikan nyata.
3. **Issue dan PR**, termasuk milik sendiri: catat bug, rencana, keterbatasan.
   Repo dengan nol issue terbaca seperti belum pernah dipakai.
4. **Pengguna di luar penulis.** Tiga rekan penulis punya ORCID; bila mereka
   memakai LAGMX dan melaporkan temuan sebagai issue, itu sinyal terkuat yang
   bisa dibangun.
5. **Arsip Zenodo** dibuat menjelang pengajuan, supaya DOI mencatat LAGMX.

### Sudah beres

editorialbot menandai empat bagian hilang di naskah — State of the field,
Software design, Research impact statement, AI usage disclosure. **Keempatnya
kini ada** di `paper.md`, bersama Validation dan Limitations.

## Langkah berikutnya

1. Buat rilis **v1.1.0** yang menandai penggantian nama. `CITATION.cff` sudah
   disiapkan di `1.1.0` / `2026-09-02`; tinggal tag dan rilisnya dibuat agar
   metadata sitasi menunjuk rilis yang benar-benar ada.

   ```bash
   gh release create v1.1.0 --title "LAGMX 1.1.0" --notes "..."
   ```
2. Jalankan poin 1–4 di atas secara bertahap sepanjang enam bulan ini.
3. Menjelang Februari 2027: buat arsip Zenodo, lalu ajukan ulang.

## Tertunda / macet

- Tag `v1.0.0`–`v1.0.3` dibuat sebelum penggantian nama dan isinya masih
  menyebut pyAutoGMX. Dibiarkan apa adanya sebagai catatan sejarah.
- Belum ada rilis dengan nama baru.

## Keputusan terakhir

- 2026-09-02 — Nama diganti **sebelum** pengajuan ulang JOSS, karena setelah
  makalah terbit nama akan terkunci di DOI dan di sitasi orang lain. Saat
  penggantian dilakukan, belum ada rujukan permanen: nihil di PyPI, nihil DOI
  arsip.
- 2026-09-02 — Tidak ada data besar dan tidak perlu Google Drive; keluaran
  simulasi dihasilkan saat dijalankan dan sudah diabaikan `.gitignore`.
- 2026-09-02 — Lisensi tetap **MIT**, berbeda dari LADOCK yang proprietary.
