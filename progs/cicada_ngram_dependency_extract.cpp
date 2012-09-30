//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

//
// dependency ngram language model learning
//

//
// input format:
// sentence1 ||| dependency1
// sentence2 ||| dependency2
// ...
//
// output format:
//
// counts/
//   root/
//     counts-xxx.gz
//     ...
//     root.idx
//   left/
//     counts-xxx.gz
//     ...
//     left.idx
//
//   left-sibling/
//     counts-xxx.gz
//     ...
//     left.idx
//
//   right/
//     counts-xxx.gz
//     ...
//     right.idx
//
//   right-sibling/
//     counts-xxx.gz
//     ...
//     right.idx
//
