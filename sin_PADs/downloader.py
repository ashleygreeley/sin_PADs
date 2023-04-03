from __future__ import annotations
from bs4 import BeautifulSoup
from typing import List
import urllib
import re
import requests
import pathlib

class Downloader:

    def __init__(self, url:str, download_dir=None) -> None:
        self.url = url
        self.download_dir=download_dir
        return

    def ls(self, match: str='*') -> List[Downloader]:

        r = requests.get(self.url)
        status_code = r.status_code
        if status_code // 100 in [4, 5]:
            raise ConnectionError(f'{self.url} returned a {status_code} error response.')

        matched_hrefs = self._search_hrefs(self.url, match=match)
        cls = type(self)
        downloaders = [None]*len(matched_hrefs)
        for i, matched_href in enumerate(matched_hrefs):
            new_url = urllib.parse.urljoin(self.url, matched_href, allow_fragments=True)
            downloaders[i] = cls(new_url, download_dir=self.download_dir)
        return downloaders

    def download(self, download_dir=None, overwrite: bool=False, stream: bool=False) -> pathlib.Path:

        if download_dir is None and self.download_dir is None:
            raise ValueError(f'download_dir kwarg needs to be set either '
                             f'in Downloader() or Downloader.download.')
        if download_dir is not None:
            self.download_dir = download_dir

        self.download_dir = pathlib.Path(self.download_dir)
        self.download_dir.mkdir(parents=True, exist_ok=True)
        file_name = pathlib.Path(self.url).name
        download_path = self.download_dir / file_name

        if (download_path.exists()) and (not overwrite):
            return download_path

        if stream:
            print('stream')
            r = requests.get(self.url, stream=True)
            file_size = int(r.headers.get('content-length'))
            downloaded_bites = 0

            megabyte = 1024 * 1024

            with open(download_path, 'wb') as f:
                for data in r.iter_content(chunk_size=5*megabyte):
                    f.write(data)
                    # Update the downloaded % in the terminal.
                    downloaded_bites += len(data)
                    download_percent = round(100 * downloaded_bites / file_size)
                    download_str = "#" * (download_percent // 5)
                    print(f'Downloading {file_name}: |{download_str:<20}| {download_percent}%', end='\r')
            print()  # Add a newline
        else:
            print('not stream')
            r = requests.get(self.url)
            with open(download_path, 'wb') as f:
                f.write(r.content)
            print(f'Downloaded {file_name}.')
        return download_path

    def name(self):
        return pathlib.Path(self.url).name

    def _search_hrefs(self, url: str, match: str='*') -> List[str]:
        request = requests.get(url)
        soup = BeautifulSoup(request.content, 'html.parser')

        match = match.replace('*', '.*')  # regex wildcard
        href_objs = soup.find_all('a', href=re.compile(match))
        # I found "?" to be in the column names so we should ignore them.
        matched_hrefs = [href_obj['href'] for href_obj in href_objs
                    if href_obj['href'].find('?')==-1]
        if len(matched_hrefs) == 0:
            raise FileNotFoundError(
                f'The url {url} does not contain any hyper '
                f'references containing the match kwarg="{match}".'
            )
        return matched_hrefs

    def __repr__(self) -> str:
        params = f'{self.url}, download_dir={self.download_dir},'
        return f'{self.__class__.__qualname__}(' + params + ')'

    def __str__(self) -> str:
        return (f'{self.__class__.__qualname__} with url={self.url} and '
               f'download_dir={self.download_dir}')
