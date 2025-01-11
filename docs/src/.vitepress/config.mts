import { defineConfig } from 'vitepress'
import { tabsMarkdownPlugin } from 'vitepress-plugin-tabs'
import mathjax3 from "markdown-it-mathjax3";
import footnote from "markdown-it-footnote";

const baseTemp = {
    base: 'REPLACE_ME_DOCUMENTER_VITEPRESS',// TODO: replace this in makedocs!
}

const navTemp = {
  nav: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
}

const nav = [
  ...navTemp.nav,
  { text: 'Tutorials', link: 'https://qutip.org/qutip-julia-tutorials/' },
  { text: 'Benchmarks', link: 'https://qutip.org/QuantumToolbox.jl/benchmarks/' },
  {
    component: 'VersionPicker'
  }
]

// https://vitepress.dev/reference/site-config
export default defineConfig({
    base: baseTemp.base,
    title: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
    description: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
    lastUpdated: true,
    cleanUrls: true,
    outDir: 'REPLACE_ME_DOCUMENTER_VITEPRESS', // This is required for MarkdownVitepress to work correctly...
    head: [
        ['link', { rel: 'icon', href: '/QuantumToolbox.jl/favicon.ico' }],
        ['link', { rel: 'icon', href: 'REPLACE_ME_DOCUMENTER_VITEPRESS_FAVICON' }],
        ['script', {src: `/QuantumToolbox.jl/versions.js`}],
        ['script', {src: `${baseTemp.base}siteinfo.js`}]
    ],
    ignoreDeadLinks: true,

    markdown: {
        math: true,

        // options for @mdit-vue/plugin-toc
        // https://github.com/mdit-vue/mdit-vue/tree/main/packages/plugin-toc#options
        toc: { level: [2, 3, 4] }, // for API page, triggered by: [[toc]]

        config(md) {
            md.use(tabsMarkdownPlugin),
                md.use(mathjax3),
                md.use(footnote)
        },
        theme: {
            light: "github-light",
            dark: "github-dark"
        }
    },
    themeConfig: {
        outline: 'deep',
        logo: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
        search: {
            provider: 'local',
            options: {
                detailedView: true
            }
        },
        nav,
        sidebar: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
        editLink: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
        socialLinks: [
            { icon: 'github', link: 'REPLACE_ME_DOCUMENTER_VITEPRESS' }
        ],
        footer: {
            message: 'Made with <a href="https://documenter.juliadocs.org/stable/" target="_blank"><strong>Documenter.jl</strong></a>, <a href="https://vitepress.dev" target="_blank"><strong>VitePress</strong></a> and <a href="https://luxdl.github.io/DocumenterVitepress.jl/stable" target="_blank"><strong>DocumenterVitepress.jl</strong></a><br>Released under the BSD 3-Clause License. Powered by the <a href="https://www.julialang.org" target="_blank">Julia Programming Language</a>.<br>',
            copyright: `© Copyright ${new Date().getUTCFullYear()} <a href="https://qutip.org/" target="_blank"><strong>QuTiP.org</strong></a>.`
        }
    }
})
